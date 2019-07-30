#!/bin/bash
source activate rdkit
gunicorn -w 4 --bind 0.0.0.0:5000 main:app
