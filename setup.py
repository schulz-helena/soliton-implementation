"""Setup script for mini-soliton-automata-software
"""

import os
from pathlib import Path

from setuptools import setup

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name = 'mini-soliton-automata',
    version = '1.0.0',
    description = 'Mini Soliton Automata Software that computes soliton paths for a given molecule and pair of exterior nodes',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author = 'Helena Schulz',
    author_email = 'schulz-helena@gmx.de',
    url = 'https://github.com/schulz-helena/soliton-implementation',
    packages = ['mini_soliton_automata', 'mini_soliton_automata.gui', 'mini_soliton_automata.res', 'mini_soliton_automata.soliton_classes', 'mini_soliton_automata.visualisations'],
    python_requires = '>3.9',
    install_requires = [
        'rdkit==2021.09.5',
        'networkx==2.7.1',
        'ffmpeg==5.0.0',
        'pysmiles==1.0.1',
        'matplotlib==3.5.1',
        'Pyqt5==5.15.4',
        'numpy==1.22.0',
        'Pillow==9.0.1'
    ],
    include_package_data=True,
    entry_points = {
        'console_scripts': [
            'mini-soliton-automata-software = mini_soliton_automata.run:main'
        ],
    },
)
