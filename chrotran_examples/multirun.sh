#!/bin/bash

./runScript.sh -t=example_1 -f=CT1di -np=16
./runScript.sh -t=example_1 -f=CT1dx -np=16
./runScript.sh -t=example_1 -f=CT1xi -np=16
./runScript.sh -t=example_1 -f=CT1xx -np=16
./runScript.sh -t=example_2 -f=CT2 -np=16

cd example_1
python CT1.py
cd ../example_2
python CT2.py

