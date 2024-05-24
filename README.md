## dssp.py
A script that will scan a .cif file and extract a sequence and its associated secondary structure. This was used over the entire [dssp database](https://pdb-redo.eu/dssp). It creates a file for each .cif file observed

Run using `python3 dssp.py source_file destination_file`
## combine.py
A quick script used to take all of the many files generated by `dssp.py` and combine them into one larger file

Run using `python3 combine.py @@@@@@@@@@@@@@@@@@@@@@@@@@@@@`
