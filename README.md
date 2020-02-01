# ChemicalStructureWebService

Web Server for Chemical Structure pictures as well as other chemical structure things

![unittest](https://github.com/mwang87/ChemicalStructureWebService/workflows/unittest/badge.svg)
![production-integration](https://github.com/mwang87/ChemicalStructureWebService/workflows/production-integration/badge.svg)

## Web API Endpoints

Full resolution of all structural information conversion

```http://localhost:5066/convert?smiles=CCO```

InChI Conversion

```http://localhost:5066/inchi?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```

InChIKey Conversion

```http://localhost:5066/inchikey?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```

Smiles Conversion

```http://localhost:5066/inchikey?smiles=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Mol Conversion

```http://localhost:5066/mol?smiles=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

ClassyFire

```http://localhost:5066/classyfire?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Formula

```http://localhost:5066/formula?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Mass

```http://localhost:5066/structuremass?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Fingerprint

```http://localhost:5066/structurefingerprint?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Image Creation

```http://localhost:5066/inchikey?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Similarity

```http://localhost:5066/structurefingerprint?inchi1=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3&smiles2=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```


## Unit Testing

``` cd test
nose2 -v
```