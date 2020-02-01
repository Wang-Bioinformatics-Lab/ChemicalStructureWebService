# Abstraction class for RDKit logic
import os
import requests
from typing import Union, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Draw import MolToFile
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprintMol

IMAGE_TYPES = ["svg", "png"]

class Molecule(object):

    def __init__(self, smiles:str=None, inchi:str=None, inchikey:str=None) -> None:
        # Baseline attributes
        # Precompute method
        self.mol = None
        self.smiles = None
        self.inchi = None
        self.inchikey = None
        # Calculate on the fly method
        # self._smiles = None
        # self._inchi = None
        # self._inchikey = None
        self._molblock = None
        self._fingerprint = None
        self._exact_mass = None
        self._formula = None

        # initial RDKit mol based on input
        # TODO: Perform fancy try catch statements here to see which one works
        if smiles:
            self._smiles = smiles
            self.mol = Chem.MolFromSmiles(smiles)
            print(self.mol)
        elif inchi:
            self._inchi = inchi
            self.mol = Chem.MolFromInchi(inchi)
        elif inchikey:
            # self._inchikey = inchikey
            self.mol = Chem.MolFromSmiles(cactus_inchikey_lookup(inchikey))

        # Get the basic properties
        self.__calc_basic_props()

    # Calculate on the fly method
    # @property
    # def smiles(self):
    #     if not self._smiles:
    #         self._smiles = Chem.MolToSmiles(self.mol)
    #     return self._smile
    
    # @property
    # def inchi(self):
    #     if not self._inchi:
    #         self._inchi = Chem.MolToInchi(self.mol)
    #     return self._inchi

    # @property
    # def inchikey(self):
    #     if not self._inchikey:
    #         self._inchikey = Chem.MolToInchiKey(self.mol)
    #     return self._inchikey


    def __calc_basic_props(self) -> None:
        """Basic properties for simple conversions"""
        try:
            self.smiles = Chem.MolToSmiles(self.mol)
            self.inchi = Chem.MolToInchi(self.mol)
            self.inchikey = Chem.MolToInchiKey(self.mol)
        except TypeError:
            pass

    def __eq__(self, other: object) -> bool:
        try:
            return self.inchikey == other.inchikey
        except AttributeError:
            return False
        except:
            raise

    def __bool__(self) -> bool:
        """Evaluate truthiness of object"""
        return bool(self.mol)

    @property
    def formula(self) -> str:
        if not self._formula:
            self._formula = CalcMolFormula(self.mol)
        return self._formula

    @property
    def molblock(self) -> str:
        if not self._molblock:
            self._molblock = Chem.MolToMolBlock(self.mol)
        return self._molblock

    @property
    def exact_mass(self) -> float:
        if not self._exact_mass:
            self._exact_mass = ExactMolWt(self.mol)
        return self._exact_mass
    
    @property
    def fingerprint(self) -> str:
        if not self._fingerprint:
            self._fingerprint = AllChem.GetMorganFingerprintAsBitVect(self.mol, radius=2, nBits=512).ToBitString()
        return self._fingerprint

    def similarity(self, mol2: 'Molecule') -> float:
        """Take another molecule and compute fingerprint similarity"""
        if not mol2 or not self:
            raise TypeError('Both molecules must be valid.')
        return FingerprintSimilarity(FingerprintMol(self.mol), FingerprintMol(mol2.mol))

    def save_image(self, path: str, height:int=250, width:int=350, imageType:str="svg") -> None:
        """Save mol SVG to given path"""
        img_type = imageType.lower()
        if img_type not in IMAGE_TYPES:
            raise TypeError(f"imageType must be in {IMAGE_TYPES}")
        if not os.path.exists(os.path.dirname(path)):
            raise OSError("Base directory does not exist")

        if img_type == "png":
            intermediate_height = 250
            intermediate_width = float(width)/float(height) * 250
            
            MolToFile(self.mol, path, size=(intermediate_width,intermediate_height), 
            subImgSize=(intermediate_width, intermediate_height), 
            fitImage=True, legends=None, imageType="svg")

            import cairosvg
            cairosvg.svg2png(url=path, write_to=path, output_width=width, output_height=height)
        elif img_type == "svg":
            MolToFile(self.mol, path, size=(width,height), 
            subImgSize=(width, height), 
            fitImage=True, legends=None, imageType="svg")


    def export_structure(self) -> dict:
        """return structure descriptors"""
        return {
            'smiles': self.smiles,
            'inchi': self.inchi,
            'inchikey': self.inchikey,
            'molblock': self.molblock
        }
        

#### HELPERS ####
# TODO: improve inchikey lookup with cactus
def cactus_inchikey_lookup(inchikey: str) -> Union[str, None]:
    url = "https://cactus.nci.nih.gov/chemical/structure/{0}/smiles".format(inchikey)
    r = requests.get(url)
    if r.ok:
        smiles = r.text.split()
        if smiles:
            return smiles[0]
    return None