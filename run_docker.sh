docker rm chemicalwebservice
docker run -d -p 5000:5000 --name chemicalwebservice chemicalwebservice /app/run_server.sh
