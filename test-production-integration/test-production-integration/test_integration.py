import requests

PRODUCTION_URL = "gnps-structure.ucsd.edu"

def test_production():
    url = f"https://{PRODUCTION_URL}/heartbeat"
    r = requests.get(url)
    r.raise_for_status()