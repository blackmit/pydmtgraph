from distutils.core import setup
from distutils.extension import Extension

PyDMTGraph = Extension(
    "PyDMTGraph",
    sources=["src/DMTGraph.cpp"],
    libraries=["boost_python39", "boost_numpy39"]
)

setup(
    name="PyDMTGraph",
    version="0.1",
    ext_modules=[PyDMTGraph]
)
