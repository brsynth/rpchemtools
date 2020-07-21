# rpchemtools

Minimalist toolbox to deal with chemicals

## Installation

```bash
conda create --name myenv python=3.6
source activate myenv
conda install --channel rdkit --channel tduigou
```

## Use
```
from chemtools.Standardizer import Standardizer
from rdkit.Chem import MolFromInchi

inchi = 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+'

params = {
    'OP_REMOVE_ISOTOPE': True,
    'OP_NEUTRALISE_CHARGE': True,
    'OP_REMOVE_STEREO': True,
    'OP_COMMUTE_INCHI': True,
    'OP_KEEP_BIGGEST': True,
    'OP_ADD_HYDROGEN': False,
    'OP_KEKULIZE': False,
    'OP_NEUTRALISE_CHARGE_LATE': True
}

mol = MolFromInchi(inchi, sanitize=False)
smol = Standardizer(sequence_fun='sequence_tunable', params=params).compute(mol)
```

## Bugs

### TD201904.01 -- stereochemistry assignment from MolToSmiles(sanitize=True) is not the same as MolToSmiles + SanitizeMol 

Source (april 2019): https://github.com/rdkit/rdkit/issues/2361

In the second case the stereo assignment is not made:

```
In [44]: m = AllChem.MolFromSmiles('[O-][n+]1onc2cc(/C=C/c3ccc(Cl)cc3)ccc21', sanitize=True)

In [45]: AllChem.MolToInchiKey(m)
Out[45]: 'AALOGNDNCMFSSI-OWOJBTEDSA-N'

In [46]: m = AllChem.MolFromSmiles('[O-][n+]1onc2cc(/C=C/c3ccc(Cl)cc3)ccc21', sanitize=False)

In [47]: AllChem.SanitizeMol(m, Chem.rdmolops.SanitizeFlags.SANITIZE_ALL)
Out[47]: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE

In [48]: AllChem.MolToInchiKey(m)
Out[48]: 'AALOGNDNCMFSSI-UHFFFAOYSA-N'
```

For the case of calling SanitizeMol after MolFromSmiles you can
force rdkit to calculate the correct InChI key by calling 
AllChem.AssignStereochemistry(m1, cleanIt=True, force=True) 
before calculating the key.
