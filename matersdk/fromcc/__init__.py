import os
from pathlib import Path


### Part . nblist
import sys
matersdk_root_dir:str = Path(__file__).parent.parent.parent.absolute()
nblist_bind_dir:str = os.path.join(matersdk_root_dir, "source", "nblist", "bind")
nblist_bind_gen_dir:str = os.path.join(nblist_bind_dir, "gen")
sys.path.append(nblist_bind_gen_dir)
# import nblist
import nblist

