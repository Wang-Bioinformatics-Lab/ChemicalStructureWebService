# ChemicalStructureWebService

Web Server for Chemical Structure pictures as well as other chemical structure things

![unittest](https://github.com/mwang87/ChemicalStructureWebService/workflows/unittest/badge.svg)
![production-integration](https://github.com/mwang87/ChemicalStructureWebService/workflows/production-integration/badge.svg)

## Web API Endpoints

Full resolution of all structural information conversion

```https://structure.gnps2.org/convert?smiles=CCO```

InChI Conversion

```https://structure.gnps2.org/inchi?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```

InChIKey Conversion

```https://structure.gnps2.org/inchikey?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```

Smiles Conversion

```https://structure.gnps2.org/inchikey?smiles=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Mol Conversion

```https://structure.gnps2.org/mol?smiles=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

ClassyFire

```https://structure.gnps2.org/classyfire?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Formula

```https://structure.gnps2.org/formula?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Mass

```https://structure.gnps2.org/structuremass?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Fingerprint

```https://structure.gnps2.org/structurefingerprint?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Image Creation

```https://structure.gnps2.org/structureimg?inchi=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3```

Structure Similarity

```https://structure.gnps2.org/structuresimilarity?inchi1=InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3&smiles2=CN1C=NC2=C1C(=O)N(C(=O)N2C)C```


## Unit Testing

``` cd test
nose2 -v
```

## Integration Testing

For local testing go to integration testing folder

```sh ./test_local.sh```

for production

```nose2 -v```
