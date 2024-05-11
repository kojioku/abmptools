# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
import subprocess
import sys

with open('README.md', 'r', encoding='utf-8') as f:
    readme = f.read()

argvs = sys.argv
if 'install' in argvs:
    try:
        subprocess.call('make')
    except:
        pass

setup(
    name='ABMPTools',
    version='1.14.0',
    description='setup tool for ABINIT-MP',
    long_description=readme,
    # install_requires=['numpy', 'pandas'],
    author='Koji Okuwaki',
    author_email='koujioku81@gmail.com',
    # license=license,
    packages=find_packages(exclude=('tests', 'docs', 'sample')),
    package_data={'abmptools': ['f90/bin/*', '../*', '../tips/*']}
)
