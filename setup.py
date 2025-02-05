import os.path
import codecs
import pathlib
from setuptools import setup, find_packages


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


# The directory containing this file
cwd = pathlib.Path(__file__).parent

# The text of the README file
README = (cwd / "README.md").read_text()

setup(
    name="py2sterimol",
    python_requires='>=3.6',
    version=get_version("py2sterimol/__init__.py"),
    description="Thin client Python interface to the original Fortran implementation of Sterimol parameters.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/JacksonBurns/py2sterimol",
    author="Jackson Burns",
    license="MIT",
    classifiers=["Programming Language :: Python :: 3"],
    install_requires=read("requirements.txt").split("\n"),
    packages=find_packages(),
    include_package_data=True,
)
