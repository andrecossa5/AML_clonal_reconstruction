# Image for variants analysis

################################################################

## Start r-base:4.2.0
FROM rocker/r-base:4.2.0

## Change as u want
LABEL mantainer = "andrea.cossa@ieo.it"

## Compilers, utilities, git
RUN apt-get --quiet update --yes \
  && apt-get --quiet upgrade --yes \
  && apt-get --quiet install --yes \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    libz-dev \
    wget \
    procps \
    git \
    cmake \
  && rm -rf /var/lib/apt/lists/*

## Conda
RUN cd / \
  && wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3 \
  && rm -f Miniconda3-latest-Linux-x86_64.sh
  
ENV PATH "$PATH:/miniconda3/bin"

## Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 
RUN tar -xf samtools-1.15.1.tar.bz2 
RUN rm samtools-1.15.1.tar.bz2 
RUN cd samtools-1.15.1 && \
    make

ENV PATH="/samtools-1.15.1/:$PATH"

## STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.5c.tar.gz
RUN tar -xzf 2.7.5c.tar.gz 
RUN rm 2.7.5c.tar.gz
RUN cd STAR-2.7.5c/source && \
    make

ENV PATH="/STAR-2.7.5c/source/STAR:$PATH"

## Conda tools and python
RUN conda install --quiet --yes -c conda-forge python=3.10 numpy pandas openjdk matplotlib seaborn
RUN pip install pysam pyclone

## Gatk
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
RUN ./gradlew bundle

## Default command to initialize the container
CMD [ "bash" ]

#################################################################
