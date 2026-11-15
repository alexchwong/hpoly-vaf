# Build environment for hpoly-vaf. Matches the toolchain the kit was authored on (Ubuntu 24.04, g++ 13).
# Build:  docker build -t hpoly-vaf-env -f kit/Dockerfile kit
# Use:    docker run --rm -it -v "$PWD":/repo -v /path/to/data:/data:ro -v /path/to/work:/work -w /repo hpoly-vaf-env
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake git ca-certificates zlib1g-dev libopenblas-dev \
      python3 python3-pip python3-pysam python3-numpy python3-scipy python3-pytest samtools tabix \
    && rm -rf /var/lib/apt/lists/*
CMD ["/bin/bash"]
