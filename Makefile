build:
	docker build -t chemicalwebservice .

server:
	docker run --rm -d -p 5065:5000 --name chemicalwebservice chemicalwebservice /app/run_server.sh

interactive:
	docker run --rm -it -p 5065:5000 --name chemicalwebservice chemicalwebservice /app/run_server.sh

dev-server:
	docker run --rm -it -p 5065:5000 --name chemicalwebservice chemicalwebservice /app/run_dev_server.sh

bash:
	docker run --rm -it -p 5065:5000 --name chemicalwebservice chemicalwebservice bash
