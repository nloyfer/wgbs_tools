# Use full image for build
FROM python:3.11-bullseye AS build

# Install poetry globally
RUN pip install poetry==1.3.2
RUN poetry --version

# Create virtualenv for deployment
ENV VIRTUAL_ENV=/opt/venv
RUN python -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install dependencies
COPY pyproject.toml poetry.lock ./
RUN poetry install --only main

WORKDIR /opt/wgbstools

COPY setup.py ./
COPY src ./src
RUN python setup.py

####################################################
FROM debian:bullseye AS samtools-build

ENV SAMTOOLS_VERSION=1.20

WORKDIR /build

RUN apt-get update && apt-get install -y wget tar bzip2 make clang zlib1g-dev libbz2-dev liblzma-dev curl
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/build && \
    make && \
    make install
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd htslib-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/build && \
    make && \
    make install

# Use slim image for deployment
FROM python:3.11-slim-bullseye

RUN apt-get update ; \
    apt-get -y install \
        samtools \
        tabix \
        bedtools \
        wget \
        unzip \
        r-base \
    ; \
    rm -rf /var/lib/apt/lists/*

WORKDIR /home/appuser

# Use the previously built virtualenv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
COPY --from=build /opt/venv /opt/venv

COPY --from=build /opt/wgbstools /opt/wgbstools
COPY --from=samtools-build /build/bin /usr/local/bin
RUN mkdir -p /opt/wgbstools/references

RUN ln -s /opt/wgbstools/src/python/wgbs_tools.py /usr/local/bin/wgbstools
RUN [ -f /opt/wgbstools/src/uxm_deconv/uxm.py ] && ln -s /opt/wgbstools/src/uxm_deconv/uxm.py /usr/local/bin/uxm || echo "uxm missing"

# App
ENTRYPOINT /bin/bash
