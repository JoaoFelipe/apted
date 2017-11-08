#!/usr/bin/env python
import os
from setuptools import setup, find_packages

try:
    import pypandoc
    LONG = pypandoc.convert('README.md', 'rst')
except (IOError, ImportError):
    LONG = "APTED algorithm for the Tree Edit Distance"


def recursive_path(pack, path):
    matches = []
    for root, _, filenames in os.walk(os.path.join(pack, path)):
        for filename in filenames:
            matches.append(os.path.join(root, filename)[len(pack) + 1:])
    return matches


setup(
    name="apted",
    version="1.0.3",
    description="APTED algorithm for the Tree Edit Distance",
    long_description=LONG,
    packages=find_packages(),
    package_data={
        "apted": (
            recursive_path("apted", "resources")
        ),
    },
    include_package_data=True,
    author=("Joao Pimentel",),
    author_email="joaofelipenp@gmail.com",
    keywords="APTED TED tree edit distance",
    url="https://github.com/JoaoFelipe/apted",
    license="MIT",
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ]
)
