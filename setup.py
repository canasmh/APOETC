#!usr/bin/env python

from setuptools import setup

setup(name='APOETC',
      version='1.0',
      packages=['APOETC'],
      scripts=['./bin/etc_script.py'],
      include_package_data=True,
      package_data={'APOETC':['../data/Sky/*.txt','../data/APO/Arctic/*.dat']}
      )

