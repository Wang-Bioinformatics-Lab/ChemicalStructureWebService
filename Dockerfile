FROM continuumio/miniconda3:latest

LABEL maintainer="Mingxun Wang mwang87@gmail.com"

WORKDIR /app
RUN apt-get update -y && apt-get install -y libxrender-dev
RUN conda create -n rdkit -c rdkit rdkit=2020.09.1.0

COPY requirements.txt /app
RUN /bin/bash -c "source activate rdkit && pip install -r requirements.txt"

COPY . /app
