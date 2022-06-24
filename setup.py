import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.0'
PACKAGE_NAME = 'syndat'
AUTHOR = 'You'
AUTHOR_EMAIL = 'nwalton1.@vols.utk.edu'
URL = 'https://github.com/Naww137/syndat'

LICENSE = 'GNU GENERAL PUBLIC LICENSE 3.0'
DESCRIPTION = 'The package creates synthetic experimental neutron cross section data.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy',
      'pandas',
      'scipy',
      'matplotlib'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )