from setuptools import setup
from ast import literal_eval


def get_version(source='xpclr/__init__.py'):
    with open(source) as sf:
        for line in sf:
            if line.startswith('__version__'):
                return literal_eval(line.split('=')[-1].lstrip())
    raise ValueError("__version__ not found")

VERSION = get_version()

DISTNAME = 'xpclr'

PACKAGE_NAME = 'xpclr'

DESCRIPTION = 'Code to compute xpclr as described in Chen 2010'

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

MAINTAINER = 'Nicholas Harding',

MAINTAINER_EMAIL = 'njh@well.ox.ac.uk',

URL = 'https://github.com/hardingnj/xpclr'

DOWNLOAD_URL = 'http://github.com/hardingnj/xpclr'

LICENSE = 'MIT'

# strictly speaking, allel requires numpy, scipy and numexpr, but numexpr
# won't install unless numpy is already installed, so leave this blank for now
# and require user to pre-install numpy, scipy and numexpr themselves
INSTALL_REQUIRES = []
CLASSIFIERS = []


def setup_package():

    metadata = dict(
        name=DISTNAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        version=VERSION,
        package_dir={'': '.'},
        packages=['xpclr'],
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
    )
    setup(**metadata)


if __name__ == '__main__':
    setup_package()