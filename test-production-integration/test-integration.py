import os
import pandas as pd
import requests
SERVER_URL = os.environ.get("SERVER_URL", "https://structure.gnps2.org")

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
        r = requests.get(url, params={"smiles" : smiles, "inchi" : inchi, "inchikey" : inchikey}, timeout=10)
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
    r = requests.get(url, params={"smiles" : smiles}, timeout=10)
    r.raise_for_status()

#def test_classyfire():
#    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
#    url = f"{SERVER_URL}/classyfire"
#    r = requests.get(url, params={"smiles" : smiles}, timeout=10)
#    r.raise_for_status()

def test_convert():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    url = f"{SERVER_URL}/convert"
    r = requests.get(url, params={"smiles" : smiles}, timeout=10)
    r.raise_for_status()
    
def test_formula():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    url = f"{SERVER_URL}/formula"
    r = requests.get(url, params={"smiles" : smiles}, timeout=10)
    r.raise_for_status()

    assert(r.text == "C8H10N4O2")

def test_input_key():
    inchikey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    url = f"{SERVER_URL}/formula"
    r = requests.get(url, params={"inchikey" : inchikey}, timeout=10)
    r.raise_for_status()

    assert(r.text == "C8H10N4O2")

def test_adduct():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    mz = 195.087
    url = f"{SERVER_URL}/adductcalc"
    r = requests.get(url, params={"smiles" : smiles, "mz": mz}, timeout=10)
    r.raise_for_status()

    assert(len(r.json()) > 5)

### Tests Relating to Formula Only API Requests ###
def test_formula_convert():
    formula = "C8H10N4O2"
    url = f"{SERVER_URL}/convert"
    r = requests.get(url, params={"formula" : formula}, timeout=10)
    
    assert(r.status_code == 400)
    assert(r.text.strip() == '{"message":"unable to import structure"}')
    
def test_formula_structuremass():
    formula = "C8H10N4O2"
    url = f"{SERVER_URL}/structuremass"
    r = requests.get(url, params={"formula" : formula}, timeout=10)
    r.raise_for_status()
    
    assert(float(r.text) >= 194.08 and float(r.text) <= 194.09)
    
def test_formula_formula():
    formula = "C8H10N4O2"
    url = f"{SERVER_URL}/formula"
    r = requests.get(url, params={"formula" : formula}, timeout=10)
    r.raise_for_status()
    
    assert(r.text.strip() == "C8H10N4O2")

def test_formula_adductcalc():
    pass
    
def test_errors_remain_handled():
    # A basic test that verifies if we put in bad strucutres, an error is still thrown
    
    # formula
    url = f"{SERVER_URL}/formula"
    
    formula = 'qwertyu'
    r = requests.get(url, params={"formula" : formula}, timeout=10)
    assert(r.status_code == 400)
    assert(r.text.strip() == '{"message":"formula cant be identified"}')
    
    smiles = 'qwertyu'
    r = requests.get(url, params={"smiles" : smiles}, timeout=10)
    assert(r.status_code == 400)
    assert(r.text.strip() == '{"message":"unable to import structure"}')
    
    # adductcalc
    url = f"{SERVER_URL}/adductcalc"
    formula = 'qwertyu'
    r = requests.get(url, params={"formula" : formula}, timeout=10)
    assert(r.status_code == 400)
    assert(r.text.strip() == '{"message":"formula cant be identified"}')
    
    smiles = 'qwertyu'
    r = requests.get(url, params={"smiles" : smiles}, timeout=10)
    assert(r.status_code == 400)
    assert(r.text.strip() == '{"message":"unable to import structure"}')
    
    # structuremass
    url = f"{SERVER_URL}/structuremass"
    formula = 'qwertyu'
    r = requests.get(url, params={"formula" : formula}, timeout=10)
    assert(r.status_code == 400)
    assert(r.text.strip() == '{"message":"formula cant be identified"}')
    
    smiles = 'qwertyu'
    r = requests.get(url, params={"smiles" : smiles}, timeout=10)
    assert(r.status_code == 400)
    assert(r.text.strip() == '{"message":"unable to import structure"}')
    