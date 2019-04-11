"""
Standardize chemicals

This is basically a rework of the standardizer.py written by Baudoin Delepine
at INRA.

@author: Baudoin Delépine, 2016-2017
@author: Thomas Duigou, 2018-2019
"""


from chemtools.Filters import Filters
from rdkit.Chem import Cleanup, SanitizeMol, SanitizeFlags


class Standardizer(object):
    """Handle standardization of compound(s) through user-defined "filters".
    """

    def __call__(self, mol):
            """Calling the Standardizer class like a function is the same
            as calling its "compute" method.
            
            Form:
                https://github.com/mcs07/MolVS/blob/master/molvs/standardize.py
            """
            return self.compute(mol)

    def __init__(self, filter_fun=None, params=None):
        """Set up parameters for the standardization
        
        :param rdmol: an RDKit Mol object
        """
        # Function to be used for standardizing compounds
        # Add you own function as method class
        if filter_fun is None:
            self.filter_fun = self.filter_minimal
        elif callable(filter_fun):  # Guess: fun_filters is the function itself
            self.filter_fun = filter_fun
        elif type(filter_fun) == str:
            self.filter_fun = globals()[filter_fun]  # Guess: fun_filters is the name of the function
        # Arguments to be passed to any custom standardization function
        self.params = params

    def filter_minimal(self, mol):
        """Minimal standardization."""
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        return mol
        
    def compute(self, mol):
        """Do the job."""
        if self.params is None:
            return self.filter_fun(mol)
        else:
            return self.filter_fun(mol, **self.params)
    
def filter_rr_legacy(mol):
    """Operations used for the first version of RetroRules
    """
    F = Filters()
    Cleanup(mol)
    SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    mol = F.remove_isotope(mol)
    mol = F.neutralise_charge(mol)
    SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    mol = F.keep_biggest(mol)
    mol = F.add_hydrogen(mol, addCoords=True)
    mol = F.kekulize(mol)
    return mol

def filter_rr_tunable(
        mol,
        OP_REMOVE_ISOTOPE=True, OP_NEUTRALISE_CHARGE=True,
        OP_REMOVE_STEREO=False, OP_COMMUTE_INCHI=False,
        OP_KEEP_BIGGEST=True, OP_ADD_HYDROGEN=True,
        OP_KEKULIZE=True
    ):
    """Tunable standardization function.
    
    Operations will made in the following order:
     1 RDKit Cleanup      -- always
     2 RDKIT SanitizeMol  -- always
     3 Remove isotope     -- optional (default: True)
     4 Neutralise charges -- optional (default: True)
     5 RDKit SanitizeMol  -- if 4 or 5
     6 Remove stereo      -- optional (default: False)
     7 Commute Inchi      -- if 6 or optional (default: False)
     8 Keep biggest       -- optional (default: True)
     9 RDKit SanitizeMol  -- if any (6, 7, 8)
    10 Add hydrogens      -- optional (default: True)
    11 Kekulize           -- optional (default: True)
    """
    F = Filters()
    # Always perform the basics..
    Cleanup(mol)
    SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    # 
    if OP_REMOVE_ISOTOPE:
        mol = F.remove_isotope(mol)
    if OP_NEUTRALISE_CHARGE:
        mol = F.neutralise_charge(mol)
    if any([OP_REMOVE_ISOTOPE, OP_REMOVE_ISOTOPE]):
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    # 
    if OP_REMOVE_STEREO:
        mol = F.remove_stereo(mol)
        OP_COMMUTE_INCHI = True
    if OP_COMMUTE_INCHI:
        mol = F.commute_inchi(mol)
    if OP_KEEP_BIGGEST:
        mol = F.keep_biggest(mol)
    if any([OP_REMOVE_STEREO, OP_COMMUTE_INCHI, OP_KEEP_BIGGEST]):
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    #
    if OP_ADD_HYDROGEN:
        mol = F.add_hydrogen(mol, addCoords=True)
    if OP_KEKULIZE:
        mol = F.kekulize(mol)
    #
    return mol
