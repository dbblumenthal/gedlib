from distutils.extension import Extension

from distutils.core import setup
from Cython.Build import cythonize


extensions = [Extension("gedlibpy",
                        sources=["gedlibpy.pyx", "src/GedLibBind.cpp"],
                        include_dirs=["include","include/lsape", "include/Eigen", "include/nomad", "include/sgtelib", "include/libsvm.3.22", "include/fann", "include/boost_1_69_0"],
                        library_dirs=["lib/fann","lib/gedlib", "lib/libsvm.3.22","lib/nomad"],
                        libraries=["doublefann","sgtelib", "svm", "nomad"],
                        language="c++",
                        extra_compile_args=["-std=c++11"],
                        extra_link_args=["-std=c++11"])]

setup(ext_modules=cythonize(extensions))


# Commande Bash : python setup.py build_ext --inplace
