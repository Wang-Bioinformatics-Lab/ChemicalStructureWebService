import requests
import math
baseurl = "http://localhost:5066/"

def unit_test(endpoint,expected_to_params):
    for (expectedVal, params) in expected_to_params:
        r = requests.get(baseurl+endpoint,params=params)
        if r.status_code == 500:
            print("error:",endpoint,"server error")
            return
        if str(r.status_code) == expectedVal:
            continue
        if expectedVal == "":
            continue
        if r.text.strip()!= expectedVal:
            print(endpoint,params,r.text,",but should get ",expectedVal)
            continue


# end point is heartbeat
def main():
    unit_test("heartbeat",[("{}",None)])
    # use ethonal
    inchi_params = {"inchi":"InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"}
    smiles_params = {"smiles":"CCO"}
    bad_inchi_params =  {"inchi":"InChI=1431H3&334213"}
    bad_smiles_params = {"smiles":"CCO**"}
    worse_smiles_params = {"smiles": "CcccC12c1"}
    smiles = "CCO"
    inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    inchikey = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    unit_test("inchikey",[(inchikey,inchi_params),
                          (inchikey,smiles_params),
                          # bad param error
                          ("400",None),
                          ("",bad_smiles_params),
                          ("None",bad_inchi_params),
                          ("400", worse_smiles_params)
                          ])
    unit_test("inchi", [(inchi,smiles_params),
                          ("400",None),
                          ("",bad_smiles_params),
                          ("400", worse_smiles_params)
                          ])

    unit_test("smiles", [(smiles,inchi_params),
                          ("400",None),
                          ("400",bad_inchi_params)])

    unit_test("mol", [    ("",inchi_params),
                          ("",smiles_params),
                          ("400",None),
                          ("400",{"smiles":"badsmiles"}),
                          ("400",{"inchi":"badinchi"})])

    unit_test("structuremass", [    ("46.041864812",inchi_params),
                              ("46.041864812",smiles_params),
                              ("400",None),
                              ("400",{"smiles":"badsmiles"}),
                              ("400",{"inchi":"badinchi"})])

    width_param = {"wdith":-100}
    height_param = {"height":-100}
    unit_test("structureimg", [    ("",inchi_params),
                              ("",smiles_params),
                              ("",smiles_params.update(width_param)),
                              ("",smiles_params.update(height_param)),
                              ("400",None),
                              ("400",{"smiles":"badsmiles"}),
                              ("400",{"inchi":"badinchi"})])

    success_test1 = {"smiles1":smiles,"inchi2":inchi}
    success_test2 = {"inchi1":inchi,"smiles2":smiles}
    success_test3 = {"smiles1":"CO","smiles2":smiles}
    fail_test1 = {"smiles1":"badsmiles","smiles2":smiles}
    fail_test2 = {"inchi1":"badinchi","smiles2":smiles}

    unit_test("structuresimilarity", [  ("1.0", success_test1),
                                  ("1.0",success_test2),
                                  ("0.3333333333333333",success_test3),
                                  ("400",fail_test1),
                                  ("400",fail_test2),
                                  ("400",None),
                                  ("400",bad_smiles_params),
                                  ("400",bad_inchi_params)])
    return


if __name__ == "__main__":
    main()
