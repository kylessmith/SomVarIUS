from setuptools import setup, Extension
from Cython.Build import cythonize
from ast import literal_eval


def versiontuple(v):
    return tuple(map(int, (v.split("."))))

try:
    import pysam
    if versiontuple(pysam.__version__) < versiontuple('0.8.1'):
        raise Exception('please upgrade pysam first, e.g.: pip install --upgrade pysam')
    from Cython.Distutils import build_ext # Cython should be installed via pysam
except ImportError:
    raise Exception('please install pysam first, e.g.: pip install --upgrade pysam')


try:
    import numpy as np
except ImportError:
    raise Exception('please install numpy first')


extensions = [Extension('SomVarIUS_calling',
                        sources=['SomVarIUS/SomVarIUS_calling.pyx'],
                        include_dirs=[np.get_include()] + pysam.get_include(),
                        define_macros=pysam.get_defines())]
                        
extensions.append(Extension('query_muts',
                        sources=['SomVarIUS/query_muts.pyx'],
                        include_dirs=[np.get_include()] + pysam.get_include(),
                        define_macros=pysam.get_defines()))


setup(
    name='SomVarIUS',
    version='1.1',
    author='Kyle S. Smith',
    author_email='kyle.s.smith@ucdenver.edu',
    url='https://github.com/kylessmith/SomVarIUS',
    description='A Python utility for calling somatic mutations from a BAM or SAM file',
    ext_modules=cythonize(extensions),
    scripts=['SomVarIUS/SomVarIUS', 'SomVarIUS/bb_fit.py','SomVarIUS/mixture_model.py', 'SomVarIUS/assign_clones.py'],
    install_requires=["pysam>=0.8.1",
                      "numpy>=1.7.0",
                      "scipy>=0.12.0",
                      "pymix"],
    classifiers=['Programming Language :: Python',
                 ],
)


