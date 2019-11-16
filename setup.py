#!usr/bin/env python

from setuptools import setup

setup(name='APOETC',
      version='0.1.0',
      author='Manuel H. Canas, Alexander Stone-Martinex, Bryson Stemock, Rogelio Ochoa, Hasan Rahman',
      author_email = 'canasmh@nmsu.edu',
      packages=['APOETC'],
      scripts=['APOETC/bin/etc_script.py'],
      include_package_data=True,
      package_data={'APOETC':['./data/APO/Arctic/*.dat','./data/APO/Filter/Johnson/*.dat','./data/APO/Filter/SDSS/*.dat']},
      description='An exposure time calculator for the ARC 3.5m telescope'
      )

