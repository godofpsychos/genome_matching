#!/bin/bash

echo "Preprocessing Started..."
python3 pre_processing.py
echo "Preprocessing Completed..."
echo "Matching Started..."
g++ processing.cpp
./a.out
echo "Matching Completed..."
echo "Results are stored in output_results.txt"
