import sys
import pandas as pd
sys.path.insert(0, "..")
import math
import Molecule


def test_molecule_parse():
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
        assert(str(expected_success) == assessment)

def test_formula():
    m = Molecule.Molecule(smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert(m.formula == "C8H10N4O2")

def test_similarity():
    m1 = Molecule.Molecule(smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    m2 = Molecule.Molecule(smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    sim = m1.similarity(m2)
    assert(sim == 1.0)

def test_image():
    m = Molecule.Molecule(smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    m.save_image(path="./caffeine_small.svg", imageType="svg")
    m.save_image(path="./caffeine_small.png", imageType="png")
    m.save_image(path="./caffeine_big.png", imageType="png", height=800, width=1200)
    m.save_image(path="./caffeine_square.png", imageType="png", height=600, width=600)