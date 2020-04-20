#    (C) Copyright 2019 Anthony D. Dutoi
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#

import os
import ctypes
import tempfile
import inspect

_local_path = os.path.dirname(__file__)



# QodeHome should be the system-local absolute path name to the Qode package.
# This will be helpful for finding include files elsewhere in the package.
QodeHome = _local_path
for _ in ("..", "..", ".."):  QodeHome = os.path.dirname(QodeHome)



class _empty(object): pass

# presently just a simple way to compile multiple files with nvcc and link them to one another 
# -- the load_C function below would have to be largely rewritten to accept C and CU files otherwise. 
# TODO: make it smarter and use caching scheme as in load_C
def load_cuda(filenames, name, flags="", include=[], directory=None, cc=None, ld=None):
    include = list(include) + [_local_path+"/PyC_types.h"]
    include = [(directory+"/"+item if item[0]!="/" else item) for item in include]
    include_dirs = {os.path.dirname(item):None for item in include}	# A cheap way to remove redundancy
    for item in include_dirs:  flags += " -I"+item
    c_files = []
    cpp_files = []
    cu_files = []
    objects=[]
    separator=" "
    for f in filenames:
        objects.append(f + ".o")
        if f.endswith(".c"):
            c_files.append(f)
        elif f.endswith(".cpp") or f.endswith(".cxx"):
            cpp_files.append(f)
        else:
            cu_files.append(f)
    if cc == None: # just compile all the files with nvcc. may not generate optimal CPU code
        cc="nvcc --shared -lcuda -lcudart --compiler-options='-std=c++11 -fPIC"+ flags + " -o "+ directory + "/" + name + " " + separator.join(filenames)
        os.system(cc)
    """
    else:
        for f in c_files:
            print("gcc -O2 -std=c99 -fPIC " + flags + " -I/usr/local/cuda/include -c -L/usr/local/cuda/lib64 -lcuda -lcudart " + " -o " + directory + "/" + f + ".o " + f)
            os.system("gcc -O2 -fPIC " + flags + " -I/usr/local/cuda/include -c -L/usr/local/cuda/lib64 -lcuda -lcudart -lstdc++ " + " -o " + directory + "/"+ f + ".o " + f)
        for f in cpp_files:
            s=""
            if f.endswith(".cpp"):
                s=".cpp"
            else:
                s=".cxx"
            print("gcc -O2 -std=c99 -fPIC -I/usr/local/cuda/include -lstdc++ -c " + flags + " -o "+ directory +"/"+ f + ".o " + f)
            os.system("gcc -O2 -std=c99 -fPIC -I/usr/local/cuda/include -lstdc++ -c " + flags + " -o "+ directory +"/"+ f + ".o " + f)
        for f in cu_files:
            print("nvcc --compiler-options='-fPIC' -c " + flags + " -o "+ directory + "/" + f + ".o " + f)
            os.system("nvcc --compiler-options='-fPIC' -c " + flags + " -o " + directory + "/"+ f + ".o " + f)
        print("gcc -O2 -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcuda -lcudart -shared -fPIC -o " + directory +"/" + name +" " + separator.join(objects))
        os.system("gcc -O2 -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcuda -lcudart -shared -fPIC -lstdc++ -o " + directory +"/" + name +" " + separator.join(objects))
    """
    return ctypes.cdll.LoadLibrary(directory+"/"+name)



def load_C(C_filestem, flags="", include=None, directory=None, cc="gcc", ld="gcc"):
    """ Loads C code found C_filestem.c, compiling if necessary into a local __Ccache__ directory """
    if directory is None:   directory = os.path.dirname(os.path.abspath(inspect.stack()[1][1]))   # assume .c code is in the same directory as the calling .py module.  See:  https://stackoverflow.com/questions/13699283/how-to-get-the-callers-filename-method-name-in-python
    if cc        is "gcc":   cc = "gcc -std=c99 -c -fPIC {} {} -o {}"
    if ld        is "gcc":    ld = "gcc -std=c99 -shared  {} {} -o {}"
    if include   is None:    include = []

    # Include the C typedefs in this directory as a dependency
    include = list(include) + [_local_path+"/"+"PyC_types.h"]	# copy include and make sure it is mutable (list not tuple)

    # Look through includes to get directories, and assume dependencies with no leading "/" are found there
    include = [(directory+"/"+item if item[0]!="/" else item) for item in include]
    include_dirs = {os.path.dirname(item):None for item in include}	# A cheap way to remove redundancy
    for item in include_dirs:  flags += " -I"+item

    # Where to put byte code
    cache = directory+"/__Ccache__/"+C_filestem
    try:     os.makedirs(cache)
    except:  pass     # A little dirty, but probably it just already exists

    # A container for associated filenames
    code = _empty()
    code.c  = directory+"/"+C_filestem+".c"
    code.h  = include	# They don't literally need to end with .h
    code.so = None

    # Look to see if there are any .so files compiled since the last source modification
    last_mod = max([os.path.getmtime(file) for file in [code.c]+code.h])
    version = 0
    for file in os.listdir(cache):
        if file.endswith("-finished.so"):
            if os.path.getmtime(cache+"/"+file) > last_mod:
                code.so = cache+"/"+file
                break	# The first one we find is good enough
            else:
                version = max(1+int(file.split("-")[0]), version)	# incrementing this when source is changed is only to aid human beings who might be cleaning up old files using shell glob syntax while newer versions are loaded and running

    # If there is no .so available, compile it, keeping in mind that other threads might be doing the same (worst case is a few extra copies)
    if code.so is None:
        version = str(version) + "-"
        tmp = tempfile.NamedTemporaryFile(dir=cache, prefix=version, delete=False)    # A thread-safe way of getting a reserved filestem
        code.o      = tmp.name+".o"
        code.tmp_so = tmp.name+".so"
        code.so     = tmp.name+"-finished.so"
        tmp.close()    # do not delete or else filestem could theoretically be reused
        print(cc.format(flags, code.c, code.o))
        os.system(cc.format(flags, code.c, code.o))
        print(ld.format(flags, code.o, code.tmp_so))
        os.system(ld.format(flags, code.o, code.tmp_so))
        os.rename(code.tmp_so, code.so)    # so that incomplete compilations from other threads will not be loaded above

    # Finally, load and return the module
    return ctypes.cdll.LoadLibrary(code.so)
