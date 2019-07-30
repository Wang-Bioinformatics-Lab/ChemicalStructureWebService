from flask import abort, jsonify, render_template, request, redirect, url_for, send_file

from app import app
import os
import random
import uuid
import string
import json
import requests

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Draw import MolToImage
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprintMol

from decorators import rdkit_handle_error

# simple test run
# input: NADA
# output: NADA
@app.route("/heartbeat")
def heartbeat():
    return "{}"

# get inchikey using either smiles or inchi
# input: smiles / inchi
# output: inchikey
@app.route("/inchikey")
@rdkit_handle_error
def inchikey():
    inchikey = ""
    if "smiles" in request.values:
        inchikey = str(Chem.MolToInchiKey(Chem.MolFromSmiles(request.values["smiles"])))
    elif "inchi" in request.values:
        inchikey = str(Chem.InchiToInchiKey(request.values["inchi"]))
    else:
        return {"message":"please input inchi or smiles"}, 400
    return inchikey

# get inchi using smiles
# input: smiles
# output: inchi
@app.route("/inchi")
@rdkit_handle_error
def inchi():
    if "smiles" in request.values:
        return str(Chem.MolToInchi(Chem.MolFromSmiles(request.values["smiles"])))
    elif "inchikey" in request.values:
        return str(Chem.MolToInchi(Chem.MolFromSmiles(cactus_inchikey_lookup(request.values["smiles"]))))
    else:
        return {"message":"please input smiles"}, 400
    

# get smiles using inchi
# input: inchi
# output: smiles
@app.route("/smiles")
@rdkit_handle_error
def smiles():
    if "inchi" not in request.values:
        return {"message":"please input inchi"}, 400
    mol = Chem.MolFromInchi(request.values["inchi"])
    if not mol:
        return {"message":"structure cant be identified"}, 400
    return str(Chem.MolToSmiles(mol))


# input: inchi / smiles
# output: mol
@app.route("/mol")
@rdkit_handle_error
def mol():
    if "smiles" in request.values:
        m = Chem.MolFromSmiles(request.values["smiles"])
    elif "inchi" in request.values:
        m = Chem.MolFromInchi(request.values["inchi"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not m:
        return {"message":"structure cant be identified"}, 400
    return Chem.MolToMolBlock(m)

# input: inchi / smiles
# output : mass (float)
# todo: keep the first token of smiles
# inchi ignore for now
@app.route("/structuremass")
@rdkit_handle_error
def structuremass():
    if "smiles" in request.values:
        m = Chem.MolFromSmiles(request.values["smiles"])
    elif "inchi" in request.values:
        m = Chem.MolFromInchi(request.values["inchi"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not m:
        return {"message":"structure cant be identified"}, 400
    return str(ExactMolWt(m))


# draw the image of structure
# rdkit
# input: smiles or inchi, width, height
@app.route("/structureimg")
@rdkit_handle_error
def structureimg():
    if "smiles" in request.values:
        m = Chem.MolFromSmiles(request.values["smiles"])
    elif "inchi" in request.values:
        m = Chem.MolFromInchi(request.values["inchi"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not m:
        return {"message":"structure cant be identified"}, 400
    #Parsing out size
    width = 350
    height = 250
    if "width" in request.values:
        try:
            if float(request.args.get('width')) > 0:
                width = float(request.args.get('width'))
        except:
            pass
    if "height" in request.values:
        try:
            if float(request.args.get('height')) > 0:
                height = float(request.args.get('height'))
        except:
            pass
    print (width,height)
    structure_images = MolToImage(m, subImgSize=(width, height), legends=None)
    output_filename = os.path.join("structure_images", str(uuid.uuid4()) + ".png")
    structure_images.save(output_filename)
    return send_file(output_filename, mimetype='image/gif')


# rdkit
@app.route("/structuresimilarity")
@rdkit_handle_error
def structuresimilarity():
    if "smiles1" in request.values:
        mol1 = Chem.MolFromSmiles(request.values["smiles1"])
    elif "inchi1" in request.values:
        mol1 = Chem.MolFromInchi(request.values["inchi1"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not mol1:
        return {"message":"unable to import structure 1."}, 400

    if "smiles2" in request.values:
        mol2 = Chem.MolFromSmiles(request.values["smiles2"])
    elif "inchi2" in request.values:
        mol2 = Chem.MolFromInchi(request.values["inchi2"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not mol2:
        return {"message":"unable to import structure 2."},400

    return str(FingerprintSimilarity(FingerprintMol(mol1),FingerprintMol(mol2)))


# ignore below functions
@app.route("/structuresimilarityjsonp")
def structuresimilarityjsonp():
    return "{}"

# from inchi, smiles, get fingerprint
# circular/Morgan fingerprint with 512 bits
@app.route("/structurefingerprint")
@rdkit_handle_error
def structurefingerprint():
    if "smiles" in request.values:
        mol = Chem.MolFromSmiles(request.values["smiles"])
    elif "inchi" in request.values:
        mol = Chem.MolFromSmiles(request.values["inchi"])
    if not mol:
        return {"message":"unable to import structure."},400
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=512).ToBitString()


def cactus_inchikey_lookup(inchikey):
    url = "https://cactus.nci.nih.gov/chemical/structure/{0}/smiles".format(inchikey)
    r = requests.get(url)
    if r.ok:
        smiles = r.text.split()
        if smiles:
            return smiles[0]
    return None

################ Old Code ####################



@app.route("/debug")
def debug():
    indigo = Indigo()
    indigo_inchi = IndigoInchi(indigo)
    m = indigo_inchi.loadMolecule("InChI=1S/C27H42O4/c1-15-7-10-27(30-14-15)16(2)24-22(31-27)12-21-19-6-5-17-11-18(28)8-9-25(17,3)20(19)13-23(29)26(21,24)4/h15-22,24,28H,5-14H2,1-4H3/t15-,16 ,17 ,18 ,19-,20 ,21 ,22 ,24 ,25 ,26-,27-/m1/s1")
    m.aromatize()

    fp = m.fingerprint("sim")

    similarity_tanimoto = indigo.similarity(fp, fp, "tanimoto")

    return str(similarity_tanimoto)

@app.route("/structuremass")
def getstructuremass():
    indigo = Indigo()

    if "inchi" in request.args:
        inchi = request.args.get('inchi')
        print("INCHI", inchi)
        inchi = inchi.lstrip()


        indigo_inchi = IndigoInchi(indigo)
        mol = indigo_inchi.loadMolecule(inchi)
    elif "smiles" in request.args:
        smiles = request.args.get('smiles')
        if "filtersalts" in request.args:
            if len(smiles.split(".")) > 1:
                if len(smiles.split(".")[0]) > smiles.split(".")[1]:
                    smiles = smiles.split(".")[0]
                else:
                    smiles = smiles.split(".")[1]

        smiles = smiles.lstrip()
        mol = indigo.loadMolecule(smiles)

    mol.dearomatize()

    return_dict = {}
    return_dict["molecularweight"] = mol.molecularWeight()
    return_dict["mostabundantmass"] = mol.mostAbundantMass()
    return_dict["monoisotopicmass"] = mol.monoisotopicMass()
    return_dict["formula"] = mol.grossFormula()

    return json.dumps(return_dict)


@app.route("/smilesstructure")
def generatesmilespng():
    indigo = Indigo()
    renderer = IndigoRenderer(indigo)

    caption = request.args.get('caption')
    if caption == None:
        caption = ""


    #Parsing out size
    width_string = request.args.get('width')
    height_string = request.args.get('height')
    width = 350
    height = 250

    try:
        if width_string != None and height_string != None:
            width = int(width_string)
            height = int(height_string)
            if width < 0:
                width = 350
            if height < 0:
                height = 250
    except ValueError:
        print("Error")


    #Determining whether input is smiles or inchi
    smiles = request.args.get('smiles')
    inchi = request.args.get('inchi')

    if (smiles == None and inchi != None):
        inchi = inchi.strip("\"")
        inchi = inchi.lstrip()
        indigo_inchi = IndigoInchi(indigo)
        mol = indigo_inchi.loadMolecule(inchi)
    else:
        smiles = smiles.strip("\"")
        smiles = smiles.lstrip();
        mol = indigo.loadMolecule(smiles)


    mol.layout() # if not called, will be done automatically by the renderer
    indigo.setOption("render-output-format", "png")
    indigo.setOption("render-comment", caption)
    indigo.setOption("render-comment-position", "top")
    indigo.setOption("render-image-size", width, height)
    indigo.setOption("render-background-color", 1.0, 1.0, 1.0)
    indigo.setOption('render-coloring', True)

    output_folder = "structure_images/"
    output_filename = str(uuid.uuid4) + ".png"

    output_path = output_folder + output_filename

    renderer.renderToFile(mol, output_path)

    return send_file(output_path, mimetype='image/gif')

@app.route("/structure_similarity/smiles")
def smilessimilarity():
    indigo = Indigo()

    structure_1 = request.args.get('structure1')
    structure_2 = request.args.get('structure2')

    m1 = indigo.loadMolecule(structure_1)
    m2 = indigo.loadMolecule(structure_2)

    #m1 = indigo.loadMolecule("CC(C)C=CCCCCC(=O)NCc1ccc(c(c1)OC)O")
    #m2 = indigo.loadMolecule("COC1=C(C=CC(=C1)C=O)O")
    # Aromatize molecules because second molecule is not in aromatic form
    m1.aromatize()
    m2.aromatize()

    fp1 = m1.fingerprint("sim");
    fp2 = m2.fingerprint("sim");

    similarity_tanimoto = indigo.similarity(fp1, fp2, "tanimoto")

    return_dict = {}
    return_dict["similarity"] = similarity_tanimoto

    return json.dumps(return_dict)

@app.route("/structure_similarity/inchi")
def inchisimilarity():
    indigo = Indigo()
    indigo_inchi = IndigoInchi(indigo);

    structure_1 = request.args.get('structure1')
    structure_2 = request.args.get('structure2')

    print(structure_1, structure_2)

    m1 = indigo_inchi.loadMolecule(structure_1)
    m2 = indigo_inchi.loadMolecule(structure_2)

    # Aromatize molecules because second molecule is not in aromatic form
    m1.aromatize()
    m2.aromatize()

    fp1 = m1.fingerprint("sim");
    fp2 = m2.fingerprint("sim");

    similarity_tanimoto = indigo.similarity(fp1, fp2, "tanimoto")

    return_dict = {}
    return_dict["similarity"] = similarity_tanimoto

    return json.dumps(return_dict)

@app.route("/inchi/inchikey")
def inchikeyinchi():
    indigo = Indigo()
    indigo_inchi = IndigoInchi(indigo);

    return indigo_inchi.getInchiKey(request.args.get('inchi'))

@app.route("/smiles/inchikey")
def inchikeysmiles():
    indigo = Indigo()
    indigo_inchi = IndigoInchi(indigo)

    smiles = request.args.get('smiles')
    m = indigo.loadMolecule(smiles)

    return indigo_inchi.getInchiKey(indigo_inchi.getInchi(m))

@app.route("/inchi/smiles")
def inchi_to_smiles():
    indigo = Indigo()
    indigo_inchi = IndigoInchi(indigo)

    inchi = request.args.get('inchi')
    inchi = inchi.lstrip()
    indigo_inchi = IndigoInchi(indigo)
    mol = indigo_inchi.loadMolecule(inchi)

    return mol.smiles()

@app.route("/smiles/inchi")
def smiles_to_inchi():
    indigo = Indigo()
    indigo_inchi = IndigoInchi(indigo)

    smiles = request.args.get('smiles')
    m = indigo.loadMolecule(smiles)

    return indigo_inchi.getInchi(m)

@app.route("/structure_similarity/smiles.jsonp")
def smilessimilarity_jsonp():
    indigo = Indigo()

    structure_1 = request.args.get('structure1')
    structure_2 = request.args.get('structure2')

    functionname = request.args.get('callback')

    m1 = indigo.loadMolecule(structure_1)
    m2 = indigo.loadMolecule(structure_2)

    #m1 = indigo.loadMolecule("CC(C)C=CCCCCC(=O)NCc1ccc(c(c1)OC)O")
    #m2 = indigo.loadMolecule("COC1=C(C=CC(=C1)C=O)O")
    # Aromatize molecules because second molecule is not in aromatic form
    m1.aromatize()
    m2.aromatize()

    fp1 = m1.fingerprint("sim");
    fp2 = m2.fingerprint("sim");

    similarity_tanimoto = indigo.similarity(fp1, fp2, "tanimoto")

    return_dict = {}
    return_dict["similarity"] = similarity_tanimoto

    return_string = functionname + "(" + str(json.dumps(return_dict)) + ");"
    return return_string



@app.route("/conversion/mol")
def convert_to_mol():
    smiles_string = request.args.get("smiles")
    m = Chem.MolFromSmiles(smiles_string)

    return Chem.MolToMolBlock(m)



if __name__ == "__main__":
    app.debug = False
    app.run(host='0.0.0.0')
