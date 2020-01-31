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

<<<<<<< HEAD
# input: inchi / smiles
# output : mass (float)
=======
# get exact mass
# input: smiles, inchi, or inchikey
# output : mass (float) as string
>>>>>>> eee3dc0c2d5cce5b820a0fbcc4c89dbba26c8542
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

# ignore below functions
@app.route("/structuresimilarityjsonp")
def structuresimilarityjsonp():
    return "{}"

# from inchi, smiles, get fingerprint
# circular/Morgan fingerprint with 512 bits
@app.route("/structurefingerprint")
@rdkit_handle_error
def structurefingerprint():
<<<<<<< HEAD
    if "smiles" in request.values:
        mol = Chem.MolFromSmiles(request.values["smiles"])
    elif "inchi" in request.values:
        mol = Chem.MolFromSmiles(request.values["inchi"])
    if not mol:
        return {"message":"unable to import structure."},400
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=512).ToBitString()

# get inchi using smiles
# input: inchikey
# output: smiles
@app.route("/inchikeyresolver")
@rdkit_handle_error
def inchikeyresolver():
    inchikey = request.values["inchikey"]
    smiles = cactus_inchikey_lookup(inchikey)

    if smiles is None:
        return "", 404

    return smiles
=======
    m = molecular_factory(request)
    if m:
        return str(m.fingerprint)
    else:
        return {"message":"unable to import structure"}, 400
>>>>>>> eee3dc0c2d5cce5b820a0fbcc4c89dbba26c8542


<<<<<<< HEAD
=======
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



>>>>>>> eee3dc0c2d5cce5b820a0fbcc4c89dbba26c8542
if __name__ == "__main__":
    app.debug = False
    app.run(host='0.0.0.0')
