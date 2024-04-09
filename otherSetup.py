from setuptools import setup

APP = ['soliton_automata/run.py']
DATA_FILES = ['soliton_automata/styles.css', 'soliton_automata/styles_m.css', 'soliton_automata/data.json']
OPTIONS = {'argv_emulation': True}

setup(
    app=APP,
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
    packages = ['soliton_automata', 'soliton_automata.gui', 'soliton_automata.res', 'soliton_automata.soliton_classes', 'soliton_automata.visualisations'],
    python_requires = '>=3.7',
    install_requires = [
        'rdkit>=2021.09.5',
        'networkx>=2.6',
        'pysmiles>=1.0.1',
        'matplotlib>=3.5.1', 
        'Pyqt5>=5.15.4',
    ],
    package_data={'': ['styles.css', 'styles_m.css', 'data.json']},
    include_package_data=True,
    entry_points = {
        'console_scripts': [
            'solitons = soliton_automata.run:main'
        ],
    }
)