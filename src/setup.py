from distutils.core import setup
from distutils.extension import Extension

dmtgraph = Extension(
    "PyDMTGraph.dmtgraph",
    sources=["PyDMTGraph/dmtgraph/DMTGraph.cpp"],
    # library_dirs=["${CONDA_PREFIX}/lib"], # uncomment if you downloaded boost with conda
    libraries=["boost_python39", "boost_numpy39"],
    extra_compile_args=["--std=c++11"],
)

setup(
    name="PyDMTGraph",
    version="0.1",
    packages=["PyDMTGraph"],
    ext_modules=[dmtgraph],
)
