from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("FittingUtilities", ["src/FittingUtilities.pyx"], include_dirs=[numpy.get_include()], extra_compile_args=["-O3", "-funroll-loops"]),],
    #package_dir = {'': 'src' }
    )
