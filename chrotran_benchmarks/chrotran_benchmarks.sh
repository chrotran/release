#!/bin/bash

cd abiotic
# ../../src/pflotran/chrotran -pflotranin abiotic.in
# python abiotic.py
# 
# cd ../abiotic_mimt
# ../../src/pflotran/chrotran -pflotranin abiotic_mimt.in
# python abiotic_mimt.py
# 
# cd ../microbe_biocide
# ../../src/pflotran/chrotran -pflotranin microbe_biocide.in
# python microbe_biocide.py
# 
# cd ../microbe_decay
# ../../src/pflotran/chrotran -pflotranin microbe_decay.in
# python microbe_decay.py
# 
# cd ../microbe_enzymatic
# ../../src/pflotran/chrotran -pflotranin microbe_enzymatic.in
# python microbe_enzymatic.py
# 
# cd ../microbe_growth
# ../../src/pflotran/chrotran -pflotranin microbe_growth.in
# python microbe_growth.py
# 
# cd ../microbe_growth_decay
# ../../src/pflotran/chrotran -pflotranin microbe_growth_decay.in
# python microbe_growth_decay.py

cd ../microbe_growth_inhibitor
../../src/pflotran/chrotran -pflotranin microbe_growth_noinhibitor.in
../../src/pflotran/chrotran -pflotranin microbe_growth_inhibitor.in
python microbe_growth_inhibitor.py

# cd ../microbe_respiration
# ../../src/pflotran/chrotran -pflotranin microbe_respiration.in
# python microbe_respiration.py
# 
# cd ..
# 
# if [ -d "results_benchmark" ]; then
#   rm -r results_benchmark
# fi
# 
# mkdir results_benchmark
# mv ./*/*.png results_benchmark
