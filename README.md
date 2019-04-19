# Toolbox to deal with chemicals

Thomas Duigou (thomas.duigou@inra.fr), INRA, 2018-2019

## Installation
```bash
pip install -e .
```

## Test
```bash
pytest
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

For the case of calling SanitizeMol after MolFromSmiles you can force rdkit to
calculate the correct InChI key by calling AllChem.AssignStereochemistry(m1, cleanIt=True, force=True)
before calculating the key.
