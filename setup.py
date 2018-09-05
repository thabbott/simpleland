from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

INC = [r'include']
LIB = [r'lib']
LINK = ['land']

cython_modules = [
    Extension('land',
              sources = ['py/land.pyx'],
              include_dirs = INC,
              library_dirs = LIB,
              libraries = LINK
              )
]

setup(ext_modules = cythonize(cython_modules))
