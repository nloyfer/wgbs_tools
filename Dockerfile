# Use full image for build
FROM python:3.11-bullseye as build

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
RUN mkdir -p /opt/wgbstools/references

RUN ln -s /opt/wgbstools/src/python/wgbs_tools.py /usr/local/bin/wgbstools
RUN [ -f /opt/wgbstools/src/uxm_deconv/uxm.py ] && ln -s /opt/wgbstools/src/uxm_deconv/uxm.py /usr/local/bin/uxm || echo "uxm missing"

# App
ENTRYPOINT /bin/bash
