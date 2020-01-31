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
from Molecule import Molecule


# This is the molecule factory that will take in a request and try to figure out what molecule it is
def molecular_factory(request) -> Molecule:
    """Given a flask request, create a molecule"""
    if "smiles" in request.values:
        m = Molecule(smiles=request.values["smiles"])
    elif "inchi" in request.values:
        m = Molecule(inchi=request.values["inchi"])
    elif "inchikey" in request.values:
        m = Molecule(inchikey=request.values["inchikey"])
    else:
        m = Molecule()
    return m


@app.route("/")
def homepage():
    return "Welcome to the chemical structure server, use our tools! This page under construction"

# simple test run
# input: NADA
# output: NADA
@app.route("/heartbeat")
def heartbeat():
    return "{}"

# Unified output in JSON format
# input: smiles, inchi, or inchikey
# output: json
@app.route("/convert")
@rdkit_handle_error
def convert():
    m = molecular_factory(request)
    if not m:
        return {"message":"unable to import structure"}, 400
    return jsonify(m.export_structure())

# get inchikey
# input: smiles, inchi, or inchikey
# output: inchikey
@app.route("/inchikey")
@rdkit_handle_error
def inchikey():
    m = molecular_factory(request)
    if m:
        return str(m.inchikey)
    else:
        return {"message":"unable to import structure"}, 400

# get classyfire using either smiles or inchi
# input: smiles / inchi
# output: classyfire
@app.route("/classyfire")
@rdkit_handle_error
def classyfire():
    inchikey = ""
    if "smiles" in request.values:
        inchikey = str(Chem.MolToInchiKey(Chem.MolFromSmiles(request.values["smiles"])))
    elif "inchi" in request.values:
        inchikey = str(Chem.InchiToInchiKey(request.values["inchi"]))
    else:
        return {"message":"please input inchi or smiles"}, 400
    
    r = requests.get("https://gnps-classyfire.ucsd.edu/entities/{}.json".format(inchikey))
    return r.text, r.status_code

# get inchi using smiles
# input: smiles
# output: inchi
@app.route("/inchi")
@rdkit_handle_error
def inchi():
    m = molecular_factory(request)
    if m:
        return str(m.inchi)
    else:
        return {"message":"unable to import structure"}, 400

# get smiles
# input: smiles, inchi, or inchikey
# output: smiles
@app.route("/smiles")
@rdkit_handle_error
def smiles():
    m = molecular_factory(request)
    if m:
        return str(m.smiles)
    else:
        return {"message":"unable to import structure"}, 400

# get molblock
# input: smiles, inchi, or inchikey
# output: molblock
@app.route("/mol")
@rdkit_handle_error
def mol():
    m = molecular_factory(request)
    if m:
        return str(m.molblock)
    else:
        return {"message":"unable to import structure"}, 400

# get exact mass
# input: smiles, inchi, or inchikey
# output : mass (float) as string
@app.route("/structuremass")
@rdkit_handle_error
def structuremass():
    m = molecular_factory(request)
    if m:
        return str(m.exact_mass)
    else:
        return {"message":"unable to import structure"}, 400

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
# input: smiles, inchi, or inchikey, width, height, imgType
@app.route("/structureimg")
@rdkit_handle_error
def structureimg():
    m = molecular_factory(request)
    if not m:
        return {"message":"unable to import structure"}, 400
    #Parsing out size
    width = int(request.args.get('width') or 350)
    height = int(request.args.get('width') or 350)
    print (width,height)
<<<<<<< HEAD
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
=======
    output_filename = os.path.join("structure_images", str(uuid.uuid4()) + ".svg")
    m.save_image(output_filename, height=height, width=width, imageType="svg")
    return send_file(output_filename, mimetype='image/svg+xml')
>>>>>>> eee3dc0c2d5cce5b820a0fbcc4c89dbba26c8542

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
        mol1 = Molecule(smiles=request.values["smiles1"])
    elif "inchi1" in request.values:
        mol1 = Molecule(inchi=request.values["inchi1"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not mol1:
        return {"message":"unable to import structure 1."}, 400

    if "smiles2" in request.values:
        mol2 = Molecule(smiles=request.values["smiles2"])
    elif "inchi2" in request.values:
        mol2 = Molecule(inchi=request.values["inchi2"])
    else:
        return {"message":"please input inchi or smiles"}, 400
    if not mol2:
        return {"message":"unable to import structure 2."},400

    return str(mol1.similarity(mol2))

# from inchi, smiles, get fingerprint
# circular/Morgan fingerprint with 512 bits
@app.route("/structurefingerprint")
@rdkit_handle_error
def structurefingerprint():
    m = molecular_factory(request)
    if m:
        return str(m.fingerprint)
    else:
        return {"message":"unable to import structure"}, 400


if __name__ == "__main__":
    app.debug = False
    app.run(host='0.0.0.0')
