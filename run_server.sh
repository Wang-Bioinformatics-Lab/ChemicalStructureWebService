#!/bin/bash

#python ./cherrypy_wrapper.py

gunicorn -w 4 --bind 0.0.0.0:5000 app:app
