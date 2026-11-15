import random, pysam, numpy as np, sys
random.seed(7); np.random.seed(7)
FO_CYCLE="TACGTACGTCTGAGCATCGATCGATGTACAGC"; NFLOW=400
FO=(FO_CYCLE*20)[:NFLOW]; KEY="TCAG"
def rnd(n):
    s=[]
    while len(s)<n:
        b=random.choice("ACGT")
        if len(s)>=1 and s[-1]==b and random.random()<0.6: continue
        s.append(b)
    return "".join(s)
left=rnd(300); right=rnd(300)
left=left[:-1]+"C"; right="G"+right[1:]
ref=left+"AAAAA"+right   # A5 at 0-based 300..305
open("ref.fa","w").write(">chrT\n"+ref+"\n"); pysam.faidx("ref.fa")
comp={"A":"T","C":"G","G":"C","T":"A"}
rc=lambda s:"".join(comp[b] for b in reversed(s))
def flowgram(seq_orient):
    t=KEY+seq_orient; i=0; fl=[]; zf=None
    for f in range(NFLOW):
        n=0
        while i<len(t) and t[i]==FO[f]:
            if i==len(KEY) and zf is None: zf=f
            i+=1; n+=1
        fl.append(n)
    assert i==len(t), "not enough flows"
    return fl,zf
hdr={"HD":{"VN":"1.4","SO":"coordinate"},"SQ":[{"SN":"chrT","LN":len(ref)}],
     "RG":[{"ID":"ABCDE.bc1","SM":"synth","PL":"IONTORRENT","FO":FO,"KS":KEY}]}
recs=[]
N=int(sys.argv[1]) if len(sys.argv)>1 else 400
for k in range(N):
    rev = (k%2==1); alt = (k%5==0)           # 20% ALT (del 1 A)
    start = 250 + (k%7) if not rev else 240 + (k%7)
    L=100
    if alt:
        seq=ref[start:300]+"AAAA"+ref[305:305+(L-(300-start)-4)]
        cig=[(0,300-start),(2,1),(0,L-(300-start))]
    else:
        seq=ref[start:start+L]; cig=[(0,L)]
    so = rc(seq) if rev else seq
    fl,zf=flowgram(so)
    zm=np.array(fl,float)*256+np.random.normal(0,18,NFLOW)
    a=pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(hdr))
    a.query_name="ABCDE:%05d:%05d"%(k+1,2*k+3); a.query_sequence=seq; a.flag=16 if rev else 0
    a.reference_id=0; a.reference_start=start; a.mapping_quality=60; a.cigar=cig
    a.query_qualities=pysam.qualitystring_to_array("5"*L)
    a.set_tag("RG","ABCDE.bc1"); a.set_tag("ZF",int(zf),"i")
    import array
    a.set_tag("ZM",array.array("h",[int(round(x)) for x in zm]))
    a.set_tag("ZP",array.array("f",[0.0,0.0,0.0]))
    recs.append((start,a,alt))
recs.sort(key=lambda r:r[0])
with pysam.AlignmentFile("reads.bam","wb",header=hdr) as o:
    for _,a,_ in recs: o.write(a)
pysam.index("reads.bam")
print("alt reads",sum(r[2] for r in recs),"of",len(recs), "ref ctx", ref[296:308])
# candidate files, written directly in TVC hotspot format
pad=ref[299]
def line(oid,r,a): return "chrT\t300\t.\t%s\t%s\t.\t.\tOID=%s;OPOS=300;OREF=%s;OALT=%s;OMAPALT=%s\n"%(r,a,oid,r,a,a)
head="##fileformat=VCFv4.1\n##allowBlockSubstitutions=true\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
open("hs1.vcf","w").write(head+line("g1.d1",pad+"A",pad))
open("hs3.vcf","w").write(head+line("g1.d1",pad+"A",pad)+line("g1.d2",pad+"AA",pad)+line("g1.i1",pad,pad+"A"))
