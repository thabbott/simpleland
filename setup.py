from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

INC = [r'include']
PYINC = [r'py/include']
LIB = [r'lib']
LINK = ['land']

cython_modules = [
    Extension('land',
              sources = ['py/src/land.pyx'],
              include_dirs = INC,
              library_dirs = LIB,
              libraries = LINK
              )
]

setup(name='CM',
      version='0.1',
      description='Python wrappers for cloud model routines',
      author='Tristan Abbott',
      author_email='tristan.h.abbott@gmail.com',
      ext_modules = cythonize(cython_modules, include_path = PYINC))
