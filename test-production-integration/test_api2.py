import requests

URL_BASE = "https://gnps-structure.ucsd.edu"

def test_smiles_to_inchikey():
    url = URL_BASE + "/inchikey?smiles=CCC(C)C1N(C)C(=O)C2CCCN2C(=O)C(OC(=O)C(Cc2ccccc2)N(C)C(=O)C(NC(=O)C(C)C(CCCCC)OC1=O)C(C)C)C(C)C"
    r = requests.get(url)
    r.raise_for_status()

    assert(r.text == "NIFSOTSGUFBSPF-UHFFFAOYSA-N")

def test_2():
    return 0