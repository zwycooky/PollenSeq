# Use a base image with the desired Linux distribution (e.g., Ubuntu)
FROM ubuntu:latest

# Set non-interactive mode for package installations
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build tools
RUN apt-get update && apt-get install -y build-essential wget nano less curl && apt-get install -y python3 && apt-get install -y wget

RUN apt-get install -y make && apt install -y gcc && apt-get install -y libz-dev && apt install -y bzip2 && apt-get install -y git

RUN apt-get install -y findutils

RUN apt-get install -y libncurses-dev && apt-get install -y liblzma-dev && apt install -y libbz2-dev && apt-get install -y parallel

RUN apt-get install -y libcurl4-openssl-dev && apt install -y unzip && apt-get install -y openjdk-17-jdk

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && tar -xzf sratoolkit.current-ubuntu64.tar.gz

RUN wget http://opengene.org/fastp/fastp && chmod a+x ./fastp

RUN git clone https://github.com/lh3/bwa.git && cd bwa && make

RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2  && bunzip2 samtools-1.17.tar.bz2 && tar -xf samtools-1.17.tar && cd samtools-1.17 && sh configure --prefix=`pwd` && make && make install

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && unzip gatk-4.4.0.0.zip && sed -i 's/python/python3/' gatk-4.4.0.0/gatk

RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar

RUN wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip

RUN rm -rf samtools-1.17.tar.bz2 gatk-4.4.0.0.zip Trimmomatic-0.39.zip

RUN apt-get install -y cpanminus

RUN cpanm -y Cwd
RUN cpanm -y Getopt::Std
RUN cpanm -y Parallel::ForkManager
RUN apt-get install -y r-base r-base-dev
# RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('HMM')"
RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN R -e "install.packages('openssl')"
RUN apt install -y  r-cran-devtools
RUN R -e "devtools::install_github('Jialab-UCR/Hapi')"
RUN R -e "install.packages('parallel')"
RUN R -e "install.packages('ASMap')"

RUN mkdir -p RephasingBins
RUN git clone https://github.com/zwycooky/PollenSeq.git
COPY ./PollenSeq/00_Softwares/Rephasing-script/ /RephasingBins/
RUN chmod +x /RephasingBins/*

#RUN git clone https://github.com/zwycooky/PollenSeq.git
# Add the R library path to the PATH variable so that it is available to your applications.
ENV PATH=$PATH:$R_HOME/lib/R/library

ENV PATH=$PATH:/RephasingBins/:/Trimmomatic-0.39/:/./:/samtools-1.17/:/gatk-4.4.0.0/:/bwa/:/sratoolkit.3.0.7-ubuntu64/bin/

# Start a new shell session when the container runs
CMD ["/bin/bash"]
