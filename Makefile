server-compose-build:
	docker-compose --compatibility build

server-compose-interactive:
	docker-compose --compatibility build
	docker-compose --compatibility --log-level DEBUG up

server-compose-server:
	docker-compose --compatibility build
	docker-compose --compatibility up -d

server-compose-production:
	docker-compose --compatibility build
	docker-compose -f docker-compose.yml -f docker-compose-production.yml --compatibility up -d

server-compose-production-interactive:
	docker-compose -f docker-compose.yml -f docker-compose-production.yml --compatibility up

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

test:
	docker exec -it chemicalwebservice_rdkit python /app/unit_test/general_test.py