from setuptools import setup, Extension
import numpy as np
from Cython.Build import cythonize
import glob

extra_compiler_args = ['-std=c++17']
depends  = [f for f in glob.glob('include/*.hpp')]
depends += [f for f in glob.glob('include/*.tpp')]
srcs     = [f for f in glob.glob('src/physics/*.cpp')]

extensions = [
    Extension(name='afterglow_ext', 
              sources=['src/afterglow_ext.pyx'] + srcs, 
              include_dirs=[np.get_include(), '.'],
              language='c++',
              depends=depends,
              extra_compile_args=extra_compiler_args)
]

setup(
    ext_modules=cythonize(extensions),
    )
