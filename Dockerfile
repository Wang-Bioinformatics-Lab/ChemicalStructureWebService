FROM continuumio/miniconda3:latest

LABEL maintainer="Mingxun Wang mwang87@gmail.com"

WORKDIR /app
RUN apt-get update -y
RUN conda create -n rdkit -c rdkit rdkit=2019.09.3.0
RUN apt-get install -y libxrender-dev
#RUN apt-get install -y inkscape

COPY requirements.txt /app
RUN /bin/bash -c "source activate rdkit && pip install -r requirements.txt"

COPY . /app
