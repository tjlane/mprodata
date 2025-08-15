#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

BIOMT_BLOCK = """REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000 -1.000000        0.00000
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file")
    args = parser.parse_args()

    input_path = Path(args.pdb_file)
    if not input_path.exists():
        raise FileNotFoundError(f"File not found: {input_path}")

    intermediate_path = input_path.with_name(f"{input_path.stem}-biomt{input_path.suffix}")

    with input_path.open("r") as f:
        lines = f.readlines()

    # Find last REMARK line index
    last_remark_idx = None
    for i, line in enumerate(lines):
        if line.startswith("REMARK"):
            last_remark_idx = i

    insert_idx = (last_remark_idx + 1) if last_remark_idx is not None else 0

    # Insert BIOMT block
    new_lines = lines[:insert_idx] + [BIOMT_BLOCK] + lines[insert_idx:]

    with intermediate_path.open("w") as f:
        f.writelines(new_lines)

    try:
        subprocess.run(
            ["phenix.pdb.biomt_reconstruction", str(intermediate_path)],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running phenix.pdb.biomt_reconstruction: {e}")
    else:
        print("phenix.pdb.biomt_reconstruction completed successfully.")

    intermediate_path.unlink()
    print(f"--> {input_path}")

if __name__ == "__main__":
    main()
