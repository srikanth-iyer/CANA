from setuptools import Extension, find_packages, setup

try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None


def readme():
    with open("README.md") as f:
        return f.read()


# Handle Cython compilation
def get_extensions():
    """
    Get the list of extensions, using Cython if available.
    """
    extensions = [
        Extension("cana.cutils", ["cana/cutils.pyx"]),
        Extension(
            "cana.canalization.cboolean_canalization",
            ["cana/canalization/cboolean_canalization.pyx"],
        ),
    ]

    if cythonize:
        # If cython is available, build from .pyx
        print("Attempting to compile from Cython.")
        return cythonize(
            extensions,
            compiler_directives={"language_level": "3"},
            include_path=[""],
        )
    else:
        # Otherwise, build from .c
        print("Cython not found. Attempting to compile from C.")
        for ext in extensions:
            ext.sources = [s.replace(".pyx", ".c") for s in ext.sources]
        return extensions


__package__ = "cana"
__description__ = "This package implements a series of methods used to study control, canalization and redundancy in Boolean networks."
__version__ = "1.0.0"

setup(
    name=__package__,
    version=__version__,
    description=__description__,
    long_description=readme(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    keywords="boolean networks canalization redundancy dynamical systems computational biology",
    url="http://github.com/rionbr/CANA",
    author="Alex Gates & Rion Brattig Correia",
    author_email="rionbr@gmail.com",
    license="MIT",
    packages=find_packages(),
    package_data={
        "datasets": [
            "cana.datasets/*.txt",
            "cana.datasets/bns/*.cnet",
            "cana.datasets/cell_collective/*.txt",
        ],
    },
    install_requires=[
        "numpy",
        "scipy",
        "networkx",
        "pandas",
        "matplotlib",
        "schematodes>=1.0.0",
    ],
    include_package_data=True,
    zip_safe=False,
    ext_modules=get_extensions(),
)