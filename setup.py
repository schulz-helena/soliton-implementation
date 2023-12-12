"""Setup script for Soliton Automata Software
"""

import os

from setuptools import setup

readme_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README.md')
with open(readme_path, "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name = 'soliton-automata-software',
    version = '2.1.0',
    description = 'Soliton Automata Software that computes soliton paths/ total legal configuration trails for a given soliton graph (and a given set of bursts)',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author = 'Helena Schulz',
    author_email = 'schulz-helena@gmx.de',
    url = 'https://github.com/schulz-helena/soliton-implementation',
    packages = ['soliton_automata', 'soliton_automata.gui', 'soliton_automata.res', 'soliton_automata.soliton_classes', 'soliton_automata.visualisations'],
    python_requires = '>=3.7',
    install_requires = [
        'rdkit>=2021.09.5',
        'networkx>=2.6',
        'pysmiles>=1.0.1',
        'matplotlib>=3.5.1',
        'Pyqt5>=5.15.4',
        'pillow==9.0.1'
    ],
    package_data={'': ['styles.css', 'styles_m.css']},
    include_package_data=True,
    entry_points = {
        'console_scripts': [
            'solitons = soliton_automata.run:main'
        ],
    },
)
