from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "sh2dwave_wrapper",
        sources=["sh2dwave_wrapper.pyx", "sh2dwave.cpp"],
        language="c++",
        include_dirs=[np.get_include()],
        extra_compile_args=["-std=c++11"],
    )
]

setup(
    name="sh2dwave_wrapper",
    ext_modules=cythonize(extensions, language_level="3"),
)