# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import subprocess
import sys

with open('README.md', 'r', encoding='utf-8') as f:
    readme = f.read()

setup(
    name='FMORmap',
    version='1.0.0',
    description='setup tool fmo calculation',
    long_description=readme,
    install_requires=[],
    author='Koji Okuwaki',
    author_email='okuwaki@rikkyo.ac.jp',
    url='http://www.cenav.org/fcews_ver1_rev2/',
    # license=license,
    packages=find_packages(exclude=('tests', 'docs', 'sample')),
)

