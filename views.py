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
from Molecule import Molecule, molecular_factory_dict
from adducts import ADDUCT_SET, get_adduct_mass

# This is the molecule factory that will take in a request and try to figure out what molecule it is
def molecular_factory(request) -> Molecule:
    """Given a flask request, create a molecule"""
    return molecular_factory_dict(request.values)

@app.route("/")
def homepage():
    return redirect('/dashinterface')


@app.route("/contributors")
def contributors():
    return render_template("contributors.html")


# simple test run
# input: NADA
# output: NADA
@app.route("/heartbeat")
def heartbeat():
    return "{}"

# Unified output in JSON format
# input: smiles, inchi, or inchikey
# output: json
@app.route("/convert", methods=['GET', 'POST'])
@rdkit_handle_error
def convert():
    m = molecular_factory(request)
    if not m.mol:
        return {"message":"unable to import structure"}, 400
    return jsonify(m.export_structure())

# get inchikey
# input: smiles, inchi, or inchikey
# output: inchikey
@app.route("/inchikey", methods=['GET', 'POST'])
@rdkit_handle_error
def inchikey():
    m = molecular_factory(request)
    if m.mol:
        return str(m.inchikey)
    else:
        return {"message":"unable to import structure"}, 400

# get classyfire using either smiles or inchi
# input: smiles / inchi
# output: classyfire
@app.route("/classyfire", methods=['GET', 'POST'])
@rdkit_handle_error
def classyfire():
    m = molecular_factory(request)
    if m:
        r = requests.get("https://classyfire.gnps2.org/entities/{}.json".format(m.inchikey))
        return r.text, r.status_code
    else:
        return {"message":"unable to import structure"}, 400

# get inchi using smiles
# input: smiles
# output: inchi
@app.route("/inchi", methods=['GET', 'POST'])
@rdkit_handle_error
def inchi():
    m = molecular_factory(request)
    if m.mol:
        return str(m.inchi)
    else:
        return {"message":"unable to import structure"}, 400

# get smiles
# input: smiles, inchi, or inchikey
# output: smiles
@app.route("/smiles", methods=['GET', 'POST'])
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
@app.route("/mol", methods=['GET', 'POST'])
@rdkit_handle_error
def mol():
    m = molecular_factory(request)
    if m.mol:
        return str(m.molblock)
    else:
        return {"message":"unable to import structure"}, 400

# get exact mass
# input: smiles, inchi, or inchikey
# output : mass (float) as string
@app.route("/structuremass", methods=['GET', 'POST'])
@rdkit_handle_error
def structuremass():
    m = molecular_factory(request)
    if m:
        return str(m.exact_mass)
    else:
        return {"message":"unable to import structure"}, 400

# input: inchi / smiles
# output : formula 
@app.route("/formula", methods=['GET', 'POST'])
@rdkit_handle_error
def formula():
    m = molecular_factory(request)
    if m:
        return str(m.formula)
    else:
        return {"message":"unable to import structure"}, 400

# input: inchi / smiles 
# output : adduct prediction 
@app.route("/adductcalc", methods=['GET', 'POST'])
@rdkit_handle_error
def calculate_adduct():
    m = molecular_factory(request)
    if m:
        exact_mass = m.exact_mass
        mz = float(request.values.get('mz'))

        all_calculations = []
        for adduct in ADDUCT_SET:
            calculated_mz, charge = get_adduct_mass(exact_mass, adduct)
            result_dict = {}
            result_dict["adduct"] = adduct
            result_dict["mz"] = calculated_mz
            result_dict["charge"] = charge
            result_dict["delta"] = mz - calculated_mz

            all_calculations.append(result_dict)
        return json.dumps(all_calculations)
    else:
        return {"message":"unable to import structure"}, 400

# draw the image of structure
# input: smiles, inchi, or inchikey, width, height, imgType
@app.route("/structureimg", methods=['GET', 'POST'])
@rdkit_handle_error
def structureimg():
    #Parsing out size
    width = int(request.args.get('width') or 350)
    height = int(request.args.get('width') or 350)
    uuid_key = str(uuid.uuid4())
    imgType = request.values.get("imgType", "png")

    m = molecular_factory(request)
    if not m.mol:
        return {"message":"unable to import structure"}, 400

    if imgType == "png":
        output_filename = os.path.join("structure_images", uuid_key + ".png")
        m.save_image(output_filename, height=height, width=width, imageType="png")
        return send_file(output_filename, mimetype='image/png')
    elif imgType == "svg":
        output_filename = os.path.join("structure_images", uuid_key + ".svg")
        m.save_image(output_filename, height=height, width=width, imageType="svg")
        return send_file(output_filename, mimetype='image/svg+xml')
    else:
        return {"message":"Please select proper file type"}, 400
    

# Calculates the structural similarity
@app.route("/structuresimilarity", methods=['GET', 'POST'])
@rdkit_handle_error
def structuresimilarity():
    smiles1 = request.values.get("smiles1", None)
    inchi1 = request.values.get("inchi1", None)

    smiles2 = request.values.get("smiles2", None)
    inchi2 = request.values.get("inchi2", None)

    mol1 = Molecule(smiles=smiles1, inchi=inchi1)
    if not mol1.mol:
        return {"message":"unable to import structure 1."}, 400

    mol2 = Molecule(smiles=smiles2, inchi=inchi2)
    if not mol2.mol:
        return {"message":"unable to import structure 2."},400

    return str(mol1.similarity(mol2))

# from inchi, smiles, get fingerprint
# circular/Morgan fingerprint with 512 bits
@app.route("/structurefingerprint", methods=['GET', 'POST'])
@rdkit_handle_error
def structurefingerprint():
    m = molecular_factory(request)
    if m.mol:
        return str(m.fingerprint)
    else:
        return {"message":"unable to import structure"}, 400

