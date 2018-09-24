from distutils.core import setup, Extension
from Cython.Build import cythonize

import numpy as np
import os
import sys


if ( str(os.environ.get('NOMAD_HOME')) == '' or str(os.environ.get('NOMAD_HOME')) == 'None'):
    print "A NOMAD_HOME environment variable is needed for building Nomad for Python (PyNomad) \n"
    exit() 

compile_args = ['-w']
link_args = []


# Look for librairies in Nomad distribution
if sys.platform.startswith('linux'):
    link_args.append('-Wl,-rpath,'+str(os.environ.get('NOMAD_HOME'))+'/lib')

# Prevent error message when changing the location of libnomad.so (OSX)
if sys.platform == 'darwin':
	link_args.append('-headerpad_max_install_names')
    
setup(
	ext_modules = cythonize(Extension(
           "PyNomad",                                # the extesion name
           sources=["PyNomad.pyx", "nomadCySimpleInterface.cpp"], # the Cython source and
	       include_dirs = [str(os.environ.get('NOMAD_HOME'))+'/src',str(os.environ.get('NOMAD_HOME'))+'/ext/sgtelib/src', np.get_include()],
           extra_compile_args=compile_args,
           extra_link_args=link_args,
           language = 'c++',
           libraries = ['nomad'],
           library_dirs = [ str(os.environ.get('NOMAD_HOME')) + '/lib'],))
)

if sys.platform == 'darwin':
      os.system('install_name_tool -change libnomad.so '+os.environ.get('NOMAD_HOME')+'/lib/libnomad.so PyNomad.so') 
