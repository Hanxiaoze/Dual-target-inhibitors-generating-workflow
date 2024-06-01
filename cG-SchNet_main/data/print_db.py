from ase import Atoms
import pickle
from ase.db import connect
import sys
import numpy as np
# 更改输出的最大宽度
#sys.setrecursionlimit(10000)
#np.set_printoptions(threshold=np.inf, linewidth=200)


# source_db="./qm9gen.db"
source_db="./world_small_30/qm9gen.db"
#source_db="./scaffold3D.db"

id = 1
with connect(source_db) as dbs:
    print(f'\n\n\nThere are 【{dbs.count()}】 molecules in this .db file !!!\n\n\n')
    for i in range(dbs.count()):
        row = dbs.get(i + 1)
        data = row.data
        print(f'\ndata {id} :\n', data)
        at = row.toatoms()
        pos = at.positions
        numbers = at.numbers
        print(f'position {id} :\n', pos)
        print(f'numbers {id} :\n', numbers)
        print(f'at {id} :\n', at)
        id += 1
    # row = dbs.get(32)
    # data = row.data
    # print(f'\ndata 32 :\n', data)
    # at = row.toatoms()
    # pos = at.positions
    # numbers = at.numbers
    # print(f'position 32 :\n', pos)
    # print(f'numbers 32 :\n', numbers)
    # print(f'at 32 :\n', at)
    print(f'\n\n\nThere are 【{dbs.count()}】 molecules in this .db file !!!\n\n\n')
