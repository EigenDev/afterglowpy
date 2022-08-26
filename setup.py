from setuptools import setup, Extension
import numpy as np
from Cython.Build import cythonize
import glob

depends  = [f for f in glob.glob('include/*.hpp')]
depends += [f for f in glob.glob('include/*.tpp')]
jet_srcs = [f for f in glob.glob('src/jet/*.cpp')]
shk_srcs = [f for f in glob.glob('src/shock/*.cpp')]
extra_compiler_args = ['-std=c++17']
test_src = ['src/jet/offaxis_struct_funcs.cpp', 'src/jet/interval.cpp']
extensions = [
    Extension(name='jet_ext', 
              sources=['src/jet_ext.pyx'] + jet_srcs + shk_srcs, 
              include_dirs=[np.get_include(), '.'],
              language='c++',
              depends=depends,
              extra_compile_args=extra_compiler_args),
    Extension(name='shock_ext', 
              sources=['src/shock_ext.pyx'] + shk_srcs + jet_srcs, 
              include_dirs=[np.get_include(), '.'],
              language='c++',
              depends=depends,
              extra_compile_args=extra_compiler_args)
]

setup(
    ext_modules=cythonize(extensions),
    )
