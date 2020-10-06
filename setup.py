try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name="ora",
    packages=["ora"],
    version="0.1",
    description="Orthologous Residue Annotator",
    author="David A. Parry",
    author_email="david.parry@igmm.ed.ac.uk",
    url='https://github.com/david-a-parry/ora',
    download_url='https://github.com/david-a-parry/ora/archive/0.1.tar.gz',
    license='MIT',
    install_requires=['requests', 'biopython'],
    scripts=["bin/ora"],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
