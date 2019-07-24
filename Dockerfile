FROM mcs07/rdkit:latest

MAINTAINER Mingxun Wang "mwang87@gmail.com"

RUN apt-get update -y
RUN apt-get install -y python3-pip python3-dev build-essential
RUN apt-get install -y default-jre

RUN pip3 install flask
RUN pip3 install gunicorn
RUN pip3 install Pillow
RUN pip3 install requests
RUN pip3 install matplotlib
RUN pip3 install cairocffi

COPY . /app
WORKDIR /app
