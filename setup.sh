#!/bin/bash

echo "Installing python3 venv and requirements, and adding project's path to PYTHONPATH ... "
python3 -m venv venv
echo 'export PYTHONPATH="$PWD:$PYTHONPATH"' >> venv/bin/activate
source venv/bin/activate
pip install -r requirements.txt