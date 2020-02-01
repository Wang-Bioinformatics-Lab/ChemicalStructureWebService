import os
import pandas as pd
import requests
SERVER_URL = os.environ.get("SERVER_URL", "https://gnps-structure.ucsd.edu")

def test_cases():
    test_cases_df = pd.read_csv("test_cases.csv", sep="\t")
    test_cases_df = test_cases_df.fillna(value="0")

    for test_case in test_cases_df.to_dict(orient="records"):
        smiles = test_case["smiles"]
        inchi = test_case["inchi"]
        inchikey = test_case["inchikey"]
        expected_success = test_case["success"]

        if smiles == "0":
            smiles = ""

        if inchi == "0":
            inchi = ""

        if inchikey == "0":
            inchikey = ""
        
        url = f"{SERVER_URL}/structuremass"
        r = requests.get(url, params={"smiles" : smiles, "inchi" : inchi, "inchikey" : inchikey})
        print(r.url, r.status_code, test_case)
        if r.status_code == 200:
            assessment = "1"
        else:
            assessment = "0"
        print(test_case, assessment)
        assert(str(expected_success) == assessment)

