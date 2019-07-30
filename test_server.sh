
#Testing inchikey conversion
http://localhost:5065/inchikey?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
http://localhost:5065/inchikey?inchi=InChI=1S/C7H15NO/c1-5-4-6(2)9-7(3)8-5/h5-8H,4H2,1-3H3

#Smiles Conversion
http://localhost:5065/smiles?inchi=InChI=1S/C7H15NO/c1-5-4-6(2)9-7(3)8-5/h5-8H,4H2,1-3H3

#InChI
http://localhost:5065/inchi?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
http://localhost:5065/inchi?inchikey=ADVPTQAUNPRNPO-UHFFFAOYSA-N 

#Structural Similarity
http://localhost:5065/structuresimilarity?smiles1=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&smiles2=CN1C=NC2=C1C(=O)N(C(=O)N2C)C

#Image Creation
http://localhost:5065/structureimg?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
http://localhost:5065/structureimg?inchi=InChI=1S/C7H15NO/c1-5-4-6(2)9-7(3)8-5/h5-8H,4H2,1-3H3

http://localhost:5065/structureimg?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&width=500&height=400

#Structure Mass
http://localhost:5065/structuremass?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
http://localhost:5065/structuremass?inchi=InChI=1S/C7H15NO/c1-5-4-6(2)9-7(3)8-5/h5-8H,4H2,1-3H3

