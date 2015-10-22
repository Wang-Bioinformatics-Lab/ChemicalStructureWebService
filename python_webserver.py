from indigo import *
from indigo_renderer import *
from indigo_inchi import *
import random
import string 
import json

from flask import Flask
from flask import request
from flask import send_file
app = Flask(__name__)

@app.route("/structuremass")
def getstructuremass():
    indigo = Indigo()
    renderer = IndigoRenderer(indigo)
    
    smiles = request.args.get('smiles')
    inchi = request.args.get('inchi')
    if (smiles == None and inchi != None):
        inchi = inchi.lstrip()
        indigo_inchi = IndigoInchi(indigo)
        mol = indigo_inchi.loadMolecule(inchi)
    else:
        filter_salts = request.args.get("filtersalts")
        if filter_salts != None:
            if len(smiles.split(".")) > 1:
                if len(smiles.split(".")[0]) > smiles.split(".")[1]:
                    smiles = smiles.split(".")[0]
                else:
                    smiles = smiles.split(".")[1]
            
        smiles = smiles.lstrip();
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
    
    print request.args
    
    try:
        if width_string != None and height_string != None:
            width = int(width_string)
            height = int(height_string)
            if width < 0:
                width = 350
            if height < 0:
                height = 250
    except ValueError:
        print "oops"
        
    
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
    output_filename = randomword(10) + ".png"
    
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
    
def randomword(length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

if __name__ == "__main__":
    app.debug = True
    app.run(host='0.0.0.0')
