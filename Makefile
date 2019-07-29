build:
	docker build -t chemicalwebservice_rdkit .

server:
	docker run --rm -d -p 5066:5000 --name chemicalwebservice_rdkit chemicalwebservice_rdkit /app/run_server.sh

interactive:
	docker run --rm -it -p 5066:5000 --name chemicalwebservice_rdkit chemicalwebservice_rdkit /app/run_server.sh

dev-server:
	docker run --rm -it -p 5066:5000 --name chemicalwebservice_rdkit chemicalwebservice_rdkit /app/run_dev_server.sh

bash:
	docker run --rm -it -p 5066:5000 --name chemicalwebservice_rdkit chemicalwebservice_rdkit bash

attach:
	docker exec -it chemicalwebservice_rdkit /bin/bash