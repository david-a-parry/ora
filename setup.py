try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import re
v_file="ora/version.py"
v_line = open(v_file, "rt").read()
v_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
match = re.search(v_re, v_line, re.M)
if match:
    verstr = match.group(1)
else:
    raise RuntimeError("Unable to find version string in {}.".format(v_file))

setup(
    name="ora",
    packages=["ora"],
    version=verstr,
    description="Orthologous Residue Annotator",
    author="David A. Parry",
    author_email="david.parry@igmm.ed.ac.uk",
    url='https://github.com/david-a-parry/ora',
    download_url='https://github.com/david-a-parry/ora/archive/{}.tar.gz'.format(verstr),
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
