import os
import json
import collections
#from pymatgen.core.structure import Molecule

def load_functional_group_data():
     """Loads functional group data from json file.
     Return list that can be easily converted into a Molecule object.
     The .json file, of course, has to be under the same directory of this function"""
     with open(os.path.join(os.path.dirname(__file__),
                        "functional_groups.json")) as f:

          data = collections.defaultdict(dict)
          #print f
          fg=[]
          for row in json.load(f):
            fg.append(row['groupname'])
            elmts=[]
            coords=[]
            els = row['sites']
            for item in els:
                elmts.append(str(item['species'][0]['element']))
                coords.append(item['xyz'])
            data[row['groupname']]=[elmts,coords]
          return data
     return None
