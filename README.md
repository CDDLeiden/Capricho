# CompoundMapper
TODO: rename the package

Initially we thought about using UniChem for mapping the compounds (hence why the package name and the unichem.py) module, but with time we've realized this was better to do with SmallWorld.

While SmallWorld requires a license for using their service, there's still a [open-source package](https://github.com/matteoferla/Python_SmallWorld_API/tree/main) for interacting with it:

The package I used for interacting with SmallWorld is a modified fork of this open-source package, so you can get the gist of the workflow. However, the added value of my fork is that it makes the requests directly to our dedicated server to SmallWorld, always retrieving the full response from the query. The `Python_SmallWorld_API`, however, often hangs without being able to retrieve the full response from the server, as the open version of it tries to limit requests to the server.

⚠️ This is not a finished version yet ⚠️

There are still a few things I would like to implement, but I just need to take the time (e.g.: being able to identify duplicate mols from fingerprints). The code is mostly in place already, I just need to make some tests and add extra functionalties.

Something else I would like to add is a command-line-interface so one can download a target-specific dataset from the terminal.

## Installation:
Clone the repo and then you can install it in your envionment by:
```bash
python -m pip install -e .
```

Disclaimer; the development of this package is still ongoing...

Usage:
```python
from CompoundMapper.unichem import UniChem
from rdkit import Chem
import pandas as pd

uchem = UniChem()

smiles = "O=C(NC[C@@H](c1ccc(Cl)cc1)N1CCOCC1)c1sc(-c2ncccn2)nc1C(F)(F)F"
inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))

res = uchem.get_connectivity(inchikey, id_type="inchikey")
df = pd.DataFrame.from_dict(res)
# df = df.query("id == 1") # This is the source_id for the ChEMBL database
comparison_cols = [c for c in df.columns if c.startswith("comparison_")]

df.assign(exact_match=df[comparison_cols].apply(lambda x: x.all(), axis=1))
```

Upcoming:
- Still to add data aggregation and filtering functions for the ChEMBL data;
- Add functions to search compounds in ChEMBL by similarity;