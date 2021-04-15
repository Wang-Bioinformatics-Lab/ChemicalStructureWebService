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

def test_image():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    url = f"{SERVER_URL}/structureimg"
    r = requests.get(url, params={"smiles" : smiles})
    r.raise_for_status()

def test_classyfire():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    url = f"{SERVER_URL}/classyfire"
    r = requests.get(url, params={"smiles" : smiles})
    r.raise_for_status()

def test_convert():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    url = f"{SERVER_URL}/convert"
    r = requests.get(url, params={"smiles" : smiles})
    r.raise_for_status()
    
def test_formula():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    url = f"{SERVER_URL}/formula"
    r = requests.get(url, params={"smiles" : smiles})
    r.raise_for_status()

    assert(r.text == "C8H10N4O2")

def test_input_key():
    inchikey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    url = f"{SERVER_URL}/formula"
    r = requests.get(url, params={"inchikey" : inchikey})
    r.raise_for_status()

    assert(r.text == "C8H10N4O2")

def test_adduct():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    mz = 195.087
    url = f"{SERVER_URL}/adductcalc"
    r = requests.get(url, params={"smiles" : smiles, "mz": mz})
    r.raise_for_status()

    assert(len(r.json()) > 5)