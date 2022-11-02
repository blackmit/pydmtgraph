import os
from platform import python_version_tuple
from distutils.core import setup
from distutils.extension import Extension

# Boost has different library names
# for different version of Python.
py_version = python_version_tuple()
libraries = [f"boost_python{py_version[0]}{py_version[1]}",
             f"boost_numpy{py_version[0]}{py_version[1]}"]

dmtgraph = Extension(
    "pydmtgraph.dmtgraph",
    sources=["pydmtgraph/dmtgraph/DMTGraph.cpp"],
    # library_dirs=["${CONDA_PREFIX}/lib"], # uncomment if you downloaded boost with conda
    libraries=libraries,
    extra_compile_args=["--std=c++11"],
)

setup(
    name="pydmtgraph",
    version="0.1",
    author="Mitchell Black",
    packages=["pydmtgraph"],
    ext_modules=[dmtgraph],
)
