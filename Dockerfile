# hash:sha256:4d7bfb7044299d17b8287dbc56136ff6bf2dfdc91a771cff72f3971a447d6de3
FROM registry.codeocean.com/codeocean/r-base:4.0.0-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && add-apt-repository -y 'deb http://deb.codeocean.com/rstudio-server-bionic/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        pandoc=1.19.2.4~dfsg-1build4 \
        rstudio-server=1.2.5033 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("WriteXLS", "4.1.0")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "PharmacoGx", \
        "S4Vectors", \
        "SummarizedExperiment" \
    ))' # Original versions: 2.0.8 0.26.1 1.18.2
