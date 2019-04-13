# Docker container to run analyses, extends GATK toolkit with a few additional
# convenience tools to identically mirror the development environment. These
# may not be necessary and could be revised in a future version for parsimony.
# inspired by Dockerfiles found at:
# https://github.com/fjukstad/seq

FROM broadinstitute/gatk:4.1.1.0

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV LANG C.UTF-8

ENV JAVA_HOME  /usr/lib/jvm/java-8-openjdk-amd64

# Set up tool dir where everything new should be installed
RUN mkdir /tools
WORKDIR /tools

USER root

#### New version bcftools and samtools Install ####
WORKDIR /tools
RUN mkdir samtools
WORKDIR samtools
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.9.tar.bz2
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar --bzip2 -xvf bcftools-1.9.tar.bz2
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar --bzip2 -xvf samtools-1.9.tar.bz2

RUN cd samtools-1.9 && \
    ./configure --prefix=/tools/samtools && \
    make && \
    make install

RUN cd bcftools-1.9 && \
    ./configure --prefix=/tools/samtools && \
    make && \
    make install

RUN rm -rf bcftools-1.9.tar.bz2 htslib-1.9.tar.bz2 samtools-1.9.tar.bz2


#### BWA Install ####
WORKDIR /tools
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.zip
RUN unzip v0.7.17.zip
RUN mv bwa-0.7.17 bwa
RUN rm v0.7.17.zip
WORKDIR bwa
RUN make


#### Picard Install ####
WORKDIR /tools
RUN mkdir picard
WORKDIR picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.18.21/picard.jar


#### Trimmomatic Install ####
WORKDIR /tools
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
RUN unzip Trimmomatic-0.38
RUN mv Trimmomatic-0.38 trimmomatic
RUN rm Trimmomatic-0.38.zip


#### Cleanup with Worker User ####
RUN groupadd -g 1001 workergroup
RUN useradd -g 1001 -u 1000 -p $(openssl passwd -1 worker) worker
RUN mkdir /home/worker
RUN mkdir /home/worker/data
WORKDIR /home/worker
RUN chmod +rx /root

#### Copy script and sample data ####
COPY process_seq.py .
COPY example.csv .
COPY data data

RUN chown -R worker:workergroup /home/worker
RUN chown -R worker:workergroup /tools

USER worker

CMD python3.6 process_seq.py /home/worker/host /home/worker/samples.csv
