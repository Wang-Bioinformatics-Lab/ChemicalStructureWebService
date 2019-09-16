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
from rdkit.Chem.Draw import MolToFile
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprintMol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

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

# input: inchi / smiles
# output : formula 
@app.route("/formula")
@rdkit_handle_error
def formula():
    if "smiles" in request.values:
        m = Chem.MolFromSmiles(request.values["smiles"])
    elif "inchi" in request.values:
        m = Chem.MolFromInchi(request.values["inchi"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not m:
        return {"message":"structure cant be identified"}, 400
    return str(CalcMolFormula(m))


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
    uuid_key = str(uuid.uuid4())
    output_svg = os.path.join("structure_images", uuid_key + ".svg")

    #TODO: The following is a hack to render at a small resolution and then scale up when we raster, it makes everything look a lot nicer
    #We should obey the aspect ratio when rendering at low resolution, but thats something we gotta do

    # structure_images = MolToFile(m, output_svg, size=(int(width), int(height)),\
    #                               subImgSize=(int(width), int(height)), \
    #                               fitImage=True, legends=None, imageType="svg")
    
    structure_images = MolToFile(m, output_svg, size=(350, 250),\
                                  subImgSize=(350, 250), \
                                  fitImage=True, legends=None, imageType="svg")

    #return send_file(output_svg, mimetype='image/svg+xml')

    #Cairo Conversion
    import cairosvg
    output_png = os.path.join("structure_images", uuid_key + ".png")
    cairosvg.svg2png(url=output_svg, write_to=output_png, output_width=width, output_height=int(250/350 * width))

    return send_file(output_png, mimetype='image/png')
    
    #Inkscape Conversion
    # output_png = os.path.join("structure_images", uuid_key + ".png")
    # cmd = "inkscape -z -e {} -w {} -h {} {}".format(output_png, width, int(250/350 * width), output_svg)
    # os.system(cmd)

    return send_file(output_png, mimetype='image/png')
    

# Calculates the structural similarity
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

if __name__ == "__main__":
    app.debug = False
    app.run(host='0.0.0.0')
