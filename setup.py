
import numpy
from distutils.core import setup, Extension

setup(name = "PLX",
      version = "1.0",
      author = "Suraj Gowda",
      author_email = "surajgowda33@gmail.com",
      py_modules = ['plx'], 
      ext_modules = [Extension("plxread", ["plxread.cpp"], include_dirs=[numpy.get_include() + '/numpy/'], language='c++')])

