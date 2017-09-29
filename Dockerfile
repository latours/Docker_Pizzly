##MAINTAINER:Sara Latour (saralatour@outlook.com)

# Base image
FROM ubuntu:16.04

# Metadata
LABEL base.image="ubuntu:16.04"
LABEL version="4"
LABEL software="Biocontainers base Image"
LABEL software.version="08252016"
LABEL description="Base image for BioDocker"
LABEL website="http://biocontainers.pro"
LABEL documentation="https://github.com/BioContainers/specs/wiki"
LABEL license="https://github.com/BioContainers/containers/blob/master/LICENSE"
LABEL tags="Genomics,Proteomics,Transcriptomics,General,Metabolomics"

# Maintainer
MAINTAINER Felipe da Veiga Leprevost <felipe@leprevost.com.br>

ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list


RUN apt-get update --yes
RUN apt-get upgrade --yes
RUN apt-get install make --yes
RUN apt-get install build-essential --yes
RUN apt-get install wget --yes
RUN apt-get install git --yes 
RUN apt-get install clang --yes
RUN apt-get purge g++ --yes
RUN apt-get install g++ --yes
RUN apt-get purge cmake --yes
RUN mkdir ~/temp
RUN cd ~/temp && wget https://cmake.org/files/v3.9/cmake-3.9.1.tar.gz && tar -xvf cmake-3.9.1.tar.gz && cd cmake-3.9.1/ && ./bootstrap && make && make install
RUN apt-get install zlib1g-dev --yes
RUN apt-get upgrade make --yes

###INSTALL THE PROGRAMS
#Download Script for Batch Processing
RUN git clone https://github.com/latours/Docker_Pizzly.git /data/batch
#Download Kallisto Program 
RUN git clone https://github.com/pmelsted/pizzly.git /data/pizzly
#Download Pizzly Program
RUN git clone https://github.com/pachterlab/kallisto.git /data/kallisto
#Build Pizzly
RUN mkdir /data/pizzly/build/ && cd /data/pizzly/build && cmake .. 
RUN cd /data/pizzly/build/ && make
RUN cd /data/pizzly/build/


###Download The GTF and FASTA Files (Reference Files)
RUN mkdir /data/reference_genome/
RUN wget http://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz > /data/reference_genome/Homo_sapiens.GRCh38.90.gtf.gz
RUN wget http://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz > /data/reference_genome/Homo_sapiens.GRCh38.cdna.all.fa.gz


RUN cp /data/batch/batch_pizzly.sh /data/pizzly/batch_pizzly.sh && chmod a+x /data/pizzly/batch_pizzly.sh

ENTRYPOINT ["/data/pizzly/batch_pizzly.sh"]
