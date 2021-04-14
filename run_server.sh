#!/bin/bash
source activate rdkit
gunicorn -w 4 --threads=2 --worker-class=gthread -b 0.0.0.0:5000 --timeout 60 --max-requests 500 --max-requests-jitter 100 --graceful-timeout 60 main:app --access-logfile /app/logs/access.log

