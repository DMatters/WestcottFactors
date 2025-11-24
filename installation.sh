#!/bin/bash

# Script to remove folders created during package-building process.
# Uninstalls westcott package, then re-installs it.
# Run this script after editing source modules in the westcott package.

rm -rf build
rm -rf dist
rm -rf westcott.egg-info

rm -f *~
rm -f westcott/*~
rm -f westcott/data_capture/*~
rm -f westcott/data_spectra/*~
rm -f tests/*~
rm -rf tests/__pycache__
rm -rf westcott/__pycache__
rm -rf __pycache__
rm -rf .tox
rm -rf notebooks/.ipynb_checkpoints/

rm -rf files_installed.txt

pip uninstall westcott<<EOF
y
EOF

#python setup.py install
# To generate list of installed files, run with the option below:
python setup.py install --record files_installed.txt

# May need to remove files manually:
#xargs rm -rf < files_installed.txt
