#!/usr/bin/bash

conda env list
# conda remove --name doublets --all
conda create --name doublets python=3.6

## Installing scrublet

conda activate doublets
pip install scrublet
pip install scanpy
