from os import path, makedirs
from dolfin.cpp.io import XDMFFile

def set_dir(name1,name2=None):
    if name2==None: x = name1
    else: x = path.join(name1, name2)
    if not path.exists(x): makedirs(x)
    return x

def set_path(name1,name2):
    return path.join(name1, name2)

def file_name(name,rs,extn):
    return name + rs + extn

def rel_error(y1,y2):
    return abs(1 - y1/y2)*100.

def save_xdmf(x, x_path):
    result = XDMFFile(x_path)
    result.parameters["flush_output"] = True
    result.parameters["functions_share_mesh"] = True
    result.write(x)
    del result, x
    return
