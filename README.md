# CompoundMapper
A collection of functionalities for mapping compounds into different IDs and retrieving data

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

smiles = "O=C(NC[C@@H](c1ccc(Cl)cc1)N1CCOCC1)c1sc(-c2ncccn2)nc1C(F)(F)F"
inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))

res = uchem.get_connectivity(inchikey, id_type="inchikey", source_id=1)
df = pd.DataFrame.from_dict(res)
# df = df.query("id == 1") # This is the source_id for the ChEMBL database
comparison_cols = [c for c in df.columns if c.startswith("comparison_")]

df.assign(exact_match=df[comparison_cols].apply(lambda x: x.all(), axis=1))
```

Upcoming:
- Still to add data aggregation and filtering functions for the ChEMBL data;
- Add functions to search compounds in ChEMBL by similarity;