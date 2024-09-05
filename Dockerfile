FROM intel/hpckit:2024.2.1-0-devel-ubuntu22.04
LABEL maintainer="Jiatai Sun <jiatai.sun@student.cup.edu.cn>"

# Installs Git.
RUN apt-get update && \
    apt-get install -y git

WORKDIR /root

# Download needed libs[ParMetis Metis GKlib]
RUN git clone https://github.com/KarypisLab/ParMETIS.git
RUN git clone https://github.com/KarypisLab/METIS.git
RUN git clone https://github.com/KarypisLab/GKlib.git

RUN cd GKlib && \
    make config openmp=set cc=icx && \
    make install -j8

RUN cd METIS && \
    make config cc=icx r64=1 && \
    make install -j8

RUN cd ParMETIS && \
    make config cc=mpiicx && \
    make install -j8
