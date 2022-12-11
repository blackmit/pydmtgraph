import os
from platform import python_version_tuple
from distutils.core import setup
from distutils.extension import Extension

# Boost has different library names
# for different version of Python.
py_version = python_version_tuple()
libraries = [f"boost_python{py_version[0]}{py_version[1]}",
             f"boost_numpy{py_version[0]}{py_version[1]}"]

# Add lib and include dirs if boost is installed in a conda environment
env_path = os.environ.get("CONDA_PREFIX")
if env_path:
    env_library_paths = (env_path, os.path.join(env_path, "Library"))
    library_dirs = [os.path.join(p, "lib") for p in env_library_paths]
    include_dirs = [os.path.join(p, "include") for p in env_library_paths]
else:
    library_dirs = None
    include_dirs = None

dmtgraph = Extension(
    "pydmtgraph.dmtgraph",
    sources=["pydmtgraph/dmtgraph/DMTGraph.cpp"],
    library_dirs=library_dirs,
    include_dirs=include_dirs,
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
