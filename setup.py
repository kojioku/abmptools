# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import subprocess
import sys

with open('README.md', 'r', encoding='utf-8') as f:
    readme = f.read()

setup(
    name='ABMPTools',
    version='1.0.0',
    description='setup tool for ABINIT-MP',
    long_description=readme,
    install_requires=[],
    author='Koji Okuwaki',
    author_email='okuwaki@rikkyo.ac.jp',
    # license=license,
    packages=find_packages(exclude=('tests', 'docs', 'sample')),
)

