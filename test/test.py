import sys
import pandas as pd
sys.path.insert(0, "..")
import math

def test_molecule():
    import Molecule
    test_cases_df = pd.read_csv("test_cases.csv", sep="\t")
    test_cases_df = test_cases_df.fillna(value="0")

    for test_case in test_cases_df.to_dict(orient="records"):
        smiles = test_case["smiles"]
        inchi = test_case["inchi"]
        inchikey = test_case["inchikey"]
        expected_success = test_case["success"]

        if smiles == "0":
            smiles = None

        if inchi == "0":
            inchi = None

        if inchikey == "0":
            inchikey = None
        
       
        m = Molecule.Molecule(smiles=smiles, inchi=inchi, inchikey=inchikey)
        assessment = ""
        if m:
            assessment = "1"
        else:
            assessment = "0"
        print(test_case, assessment)
        assert(str(expected_success) == assessment)