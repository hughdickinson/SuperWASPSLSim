from setuptools import setup

# from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np
import os

basedir = os.path.split(os.getcwd())[0]

extensions = [
    Extension(
        "PyNonUniformCircularSource",
        ["PyNonUniformCircularSource.pyx"],
        include_dirs=[
            os.path.join(basedir, "include"),
            os.path.join(basedir, "include", "boost"),
        ],
        extra_compile_args=["-std=c++14"],
    )
]

setup(
    name="py-nonuniform-circular-source",
    version="0.1",
    author="Hugh Dickinson",
    author_email="hugh.dickinson@open.ac.uk",
    ext_modules=cythonize(extensions),
)
