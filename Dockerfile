FROM continuumio/miniconda3:latest

MAINTAINER Mingxun Wang "mwang87@gmail.com"

WORKDIR /app
COPY requirements.txt /app
RUN apt-get update -y
RUN conda create -n rdkit -c rdkit rdkit
RUN apt-get install -y libxrender-dev

RUN /bin/bash -c "source activate rdkit && pip install -r requirements.txt"

COPY . /app
