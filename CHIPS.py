# CHIPs — Chromatographic Hydrophobicity Index for Peptides
# Copyright (C) 2024 Oleg V. Krokhin, Alexandre Préfontaine
# University of Manitoba, Manitoba Centre for Proteomics and Systems Biology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
Program: CHIPS.py
Author: Alexandre Prefontaine
Affiliation: University of Manitoba, Manitoba Centre for Proteomics and Systems Biology
Date of Last Update: 2024-06-20

Command-line interface for CHIPs (Chromatographic Hydrophobicity Index for Peptides).

Usage:
    python CHIPS.py <input.csv> <FA|TFA>

Arguments:
    input.csv   Single-column CSV file containing peptide sequences.
                A header row is detected and skipped automatically.
    FA|TFA      Ion pairing agent used in the HPLC gradient.

Output:
    A CSV file in the same directory as the input, named:
    <input>_predicted_HI.csv

Example:
    python CHIPS.py my_peptides.csv FA
"""

import sys
import csv
from pathlib import Path
from CHIPS_core import calculate_HI, VALID_AA, KNOWN_HEADERS


def is_valid_sequence(value: str) -> bool:
    """
    Returns True if the value looks like a peptide sequence.
    Rejects known header strings and sequences containing invalid characters.
    """
    cleaned = value.upper().strip()
    if cleaned.lower() in KNOWN_HEADERS:
        return False
    return all(c in VALID_AA for c in cleaned)


def parse_args():
    """
    Reads and validates command line arguments.
    Returns a tuple of (csv_path, phase).
    """
    if len(sys.argv) != 3:
        print("Usage: python CHIPS.py <input.csv> <FA|TFA>")
        sys.exit(1)

    csv_path = sys.argv[1]
    phase    = sys.argv[2].upper().strip()

    if phase not in ('FA', 'TFA'):
        print(f"Error: phase must be FA or TFA, got '{phase}'")
        sys.exit(1)

    return csv_path, phase


def main():
    # 1. Get and validate arguments
    csv_path, phase = parse_args()

    # 2. Build output filename in the same directory as the input
    input_path  = Path(csv_path)
    output_path = input_path.parent / (input_path.stem + '_predicted_HI.csv')

    # 3. Read, compute, write
    with open(input_path, 'r', newline='') as infile, \
         open(output_path, 'w', newline='') as outfile:

        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        writer.writerow(['Peptide', 'HI'])

        for row in reader:
            if not row:
                continue
            sequence = row[0].upper().strip()
            if not is_valid_sequence(sequence):
                continue
            try:
                hi = calculate_HI(sequence, phase)
                writer.writerow([sequence, hi])
            except Exception as e:
                print(f"Skipping '{sequence}': {e}")

    print(f"Done. Results written to {output_path}")


if __name__ == "__main__":
    main()
