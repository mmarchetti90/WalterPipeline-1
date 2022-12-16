FROM continuumio/miniconda3:4.12.0

### UPDATING CONDA ------------------------- ###

RUN conda update -y conda

### INSTALLING PIPELINE PACKAGES ----------- ###

# Adding bioconda to the list of channels
RUN conda config --add channels bioconda

# Adding conda-forge to the list of channels
RUN conda config --add channels conda-forge

# Installing mamba
RUN conda install -y mamba

# Installing packages
# N.B. have to add weasyprint and specify tb-profiler=4.4.0 or the version automatically installed is 2.8
RUN mamba install -y \
    bbmap \
    bcftools \
    bowtie2 \
    bwa \
    cutadapt \
    fastqc \
    gatk \
    kraken2 \
    lofreq \
    mummer \
    mykrobe \
    picard \
    sambamba \
    tabix \
    tb-profiler=4.4.0 \
    trim-galore \
    weasyprint && \
    conda clean -afty

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]