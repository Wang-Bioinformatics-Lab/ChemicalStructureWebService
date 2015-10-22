from cherrypy import wsgiserver
from python_webserver import app

d = wsgiserver.WSGIPathInfoDispatcher({'/': app})
server = wsgiserver.CherryPyWSGIServer(('0.0.0.0', 5000), d)

if __name__ == '__main__':
   try:
      server.start()
   except KeyboardInterrupt:
      server.stop()
