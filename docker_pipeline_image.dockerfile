FROM continuumio/miniconda3:4.9.2

### UPDATING CONDA ------------------------- ###

RUN conda update -y conda

### INSTALLING PIPELINE PACKAGES ----------- ###

# Adding bioconda to the list of channels
RUN conda config --add channels bioconda

# Installing mamba
RUN conda install -c conda-forge -y mamba

RUN mamba install -c conda-forge -y \
    bbmap \
    bcftools \
    bowtie2 \
    bwa \
    cutadapt \
    fastqc \
    gatk \
    kraken \
    lofreq \
    mummer \
    mykrobe \
    picard \
    sambamba \
    tabix \
    tb-profiler \
    trim-galore && \
    conda clean -afty

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]