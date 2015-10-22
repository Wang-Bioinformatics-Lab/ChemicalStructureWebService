#!/usr/bin/python


import sys
import getopt
from indigo import *
from indigo_renderer import *
from indigo_inchi import *

def usage():
    print "<input mgf>"
    

def smiles_renderer():
    indigo = Indigo()
    renderer = IndigoRenderer(indigo)
    
    mol = indigo.loadMolecule("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    mol.layout() # if not called, will be done automatically by the renderer 
    indigo.setOption("render-output-format", "png")
    indigo.setOption("render-comment", "Caffeine")
    indigo.setOption("render-comment-position", "top")
    indigo.setOption("render-image-size", 200, 250)
    indigo.setOption("render-background-color", 1.0, 1.0, 1.0)
    renderer.renderToFile(mol, "caffeine.png")
    
def inchi_renderer():
    indigo = Indigo()
    indigo_inchi = IndigoInchi(indigo)
    renderer = IndigoRenderer(indigo)
    
    mol = indigo_inchi.loadMolecule("InChI=1S/C7H15NO/c1-5-4-6(2)9-7(3)8-5/h5-8H,4H2,1-3H3")
    mol.layout() # if not called, will be done automatically by the renderer 
    indigo.setOption("render-output-format", "png")
    indigo.setOption("render-comment", "Caffeine")
    indigo.setOption("render-comment-position", "top")
    indigo.setOption("render-image-size", 200, 250)
    indigo.setOption("render-background-color", 1.0, 1.0, 1.0)
    renderer.renderToFile(mol, "inchi.png")

    
def main():
    inchi_renderer()
        
    
    
if __name__ == "__main__":
    main()
