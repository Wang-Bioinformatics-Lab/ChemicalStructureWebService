"""Test the speed of server with different backend designs

1) mwang = 18.703308335971087,18.248935278039426,18.457531138090417 seconds to execute 500 calls
2) JvS - precalculated = 21.32231805101037,20.969998513930477,20.78585294005461 seconds to execute 500 calls
3) JvS - calc on fly (+cache) = 18.714338005986065,18.474315379047766,18.3609054380795 seconds to execute 500 calls
"""
import requests
import time
baseurl = "http://localhost:5066/"

MAX_CALLS = 500

def main():
    start = time.perf_counter()
    for x in range(MAX_CALLS):
        r = requests.get(baseurl+'inchikey', params={'smiles': 'C'*x})
        print(r.text.strip())
    end = time.perf_counter()
    print(f"Took {end-start} seconds to execute {MAX_CALLS} calls")

if __name__ == '__main__':
    main()