FROM ubuntu:latest
MAINTAINER Mingxun Wang "mwang87@gmail.com"

RUN apt-get update -y
RUN apt-get install -y python-pip python-dev build-essential
RUN apt-get install -y default-jre
RUN pip install flask
RUN pip install gunicorn


COPY . /app
WORKDIR /app
