cmake_minimum_required(VERSION 3.16)
project(tvc_standalone C CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)
set(TS "" CACHE PATH "TS checkout")
set(DEPS "" CACHE PATH "deps dir")
set(LAPACK_LINK "openblas" CACHE STRING "LAPACK/BLAS link items")
set(EXTRA_SRC "" CACHE STRING "extra sources (sandbox LAPACK shim only)")
set(A ${TS}/Analysis)
set(V ${A}/VariantCaller)
set(X ${TS}/external)
if(DEFINED ENV{CONDA_PREFIX})
  include_directories($ENV{CONDA_PREFIX}/include)
  link_directories($ENV{CONDA_PREFIX}/lib)
endif()
add_definitions(-DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -DALIGNSTATS_IGNORE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse3 -fpermissive -w -include unistd.h -include cmath -include cstdint -idirafter ${A}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -w")
include_directories(${CMAKE_SOURCE_DIR}/../gen ${DEPS}/bt/include/bamtools ${DEPS}/RcppArmadillo/inst/include ${DEPS}/bh/inst/include
  ${X}/vcflib ${X}/freebayes/src ${X}/jsoncpp-src-amalgated0.6.0-rc1
  ${V} ${V}/EnsembleEval ${V}/Bookkeeping ${V}/Reads ${V}/Splice ${V}/Filter ${V}/Consensus ${V}/tvcutils ${V}/IndelAssembly
  ${A}/realignment ${A}/Calibration ${A}/BaseCaller ${A}/Util ${TS}/Analysis/Image ${A}/file-io ${A}/AnalysisOrg ${A}/AnalysisOrg/IO)
set(COMMON
  ${V}/BAMWalkerEngine.cpp ${V}/OrderedBAMWriter.cpp ${V}/SampleManager.cpp ${V}/TargetsManager.cpp ${V}/HandleVariant.cpp
  ${V}/HotspotReader.cpp ${V}/MetricsManager.cpp ${V}/HypothesisEvaluator.cpp ${V}/MolecularTag.cpp
  ${V}/Bookkeeping/MiscUtil.cpp ${V}/Bookkeeping/ExtendParameters.cpp ${V}/Bookkeeping/InputStructures.cpp ${V}/Bookkeeping/VcfFormat.cpp
  ${V}/Reads/ExtendedReadInfo.cpp
  ${V}/Splice/ErrorMotifs.cpp ${V}/Splice/LocalContext.cpp ${V}/Splice/ClassifyVariant.cpp ${V}/Splice/SpliceVariantHypotheses.cpp
  ${V}/Filter/DecisionTreeData.cpp ${V}/Filter/VariantAssist.cpp
  ${V}/EnsembleEval/BiasGenerator.cpp ${V}/EnsembleEval/SigmaGenerator.cpp ${V}/EnsembleEval/SkewGenerator.cpp
  ${V}/EnsembleEval/PosteriorInference.cpp ${V}/EnsembleEval/ShortStack.cpp ${V}/EnsembleEval/StackEngine.cpp ${V}/EnsembleEval/CrossHypotheses.cpp
  
  ${X}/vcflib/Variant.cpp ${X}/vcflib/split.cpp ${X}/vcflib/tabixpp/tabix.cpp ${X}/vcflib/tabixpp/index.c ${X}/vcflib/tabixpp/bgzf.c
  ${X}/vcflib/smithwaterman/LeftAlign.cpp ${X}/vcflib/smithwaterman/Repeats.cpp ${X}/vcflib/smithwaterman/IndelAllele.cpp ${X}/vcflib/smithwaterman/SmithWatermanGotoh.cpp
  ${X}/freebayes/src/AlleleParser.cpp
  ${A}/BaseCaller/PIDloop.cpp ${A}/BaseCaller/DPTreephaser.cpp ${A}/BaseCaller/MolecularTagTrimmer.cpp ${A}/BaseCaller/TreephaserSSE.cpp
  ${A}/Calibration/LinearCalibrationModel.cpp ${A}/Calibration/FlowAlignment.cpp
  ${A}/Util/OptArgs.cpp ${A}/AnalysisOrg/IO/OptBase.cpp ${A}/realignment/Realigner.cpp
  ${X}/jsoncpp-src-amalgated0.6.0-rc1/jsoncpp.cpp ${CMAKE_SOURCE_DIR}/../gen/IonVersion.cpp ${A}/file-io/ion_util.c ${A}/file-io/ion_error.c ${A}/file-io/ion_alloc.c ${EXTRA_SRC})
add_executable(tvc ${V}/VariantCaller.cpp ${COMMON})
target_link_libraries(tvc ${DEPS}/bt/lib/libbamtools.a z pthread ${LAPACK_LINK})
