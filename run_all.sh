#!/bin/bash
python3 extend.py
python3 scan_pfam.py blastp_result.fasta hmmscan_result_big.csv
python3 scan_pfam.py input-z2.fasta hmmscan_result_small.csv
python3 fisher.py hmmscan_result_big.csv hmmscan_result_small.csv
