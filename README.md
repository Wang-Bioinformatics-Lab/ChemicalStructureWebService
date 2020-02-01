# ChemicalStructureWebService

Web Server for Chemical Structure pictures as well as other chemical structure things

![unittest](https://github.com/mwang87/ChemicalStructureWebService/workflows/unittest/badge.svg)
![production-integration](https://github.com/mwang87/ChemicalStructureWebService/workflows/production-integration/badge.svg)

## Web API Endpoints

Full resolution of all structural information conversion

```https://gnps-structure.ucsd.edu/convert?smiles=CCO```

InChI Conversion

```https://gnps-structure.ucsd.edu/inchi?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```

InChIKey Conversion

```https://gnps-structure.ucsd.edu/inchikey?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```

Smiles Conversion

```https://gnps-structure.ucsd.edu/inchikey?smiles=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Mol Conversion

```https://gnps-structure.ucsd.edu/mol?smiles=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

ClassyFire

```https://gnps-structure.ucsd.edu/classyfire?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Formula

```https://gnps-structure.ucsd.edu/formula?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Mass

```https://gnps-structure.ucsd.edu/structuremass?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Fingerprint

```https://gnps-structure.ucsd.edu/structurefingerprint?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Image Creation

```https://gnps-structure.ucsd.edu/structureimg?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Similarity

```https://gnps-structure.ucsd.edu/structuresimilarity?inchi1=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3&smiles2=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```


## Unit Testing

``` cd test
nose2 -v
```

## Integration Testing

For local testing go to integration testing folder

sh ./run_local.sh

for production

nose2 -v