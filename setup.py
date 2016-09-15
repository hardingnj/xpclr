from setuptools import setup

setup(name='xpclr',
      version='0.8.0',
      description='Code to compute xpclr as described in Chen 2010',
      url='http://github.com/hardingnj/xpclr',
      author='Nicholas Harding',
      author_email='njh@well.ox.ac.uk',
      license='MIT',
      packages=['xpclr'],
      scripts=["bin/compute_xpclr.py"],
      zip_safe=False)