#!/usr/bin/env python3
"""
Generate mock FASTQ files for testing.
Creates minimal synthetic single-end FASTQ data.
"""

import gzip
import random
import sys
from pathlib import Path

def make_mock_fastq(outfile: Path, n_reads: int = 1000,
                    read_length: int = 75, seed: int = 42):
    """Generate a synthetic FASTQ.gz with random reads."""
    random.seed(seed)
    bases = "ACGT"
    quals = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"

    with gzip.open(outfile, "wt") as f:
        for i in range(n_reads):
            seq  = "".join(random.choices(bases, k=read_length))
            qual = quals[:read_length]
            f.write(f"@read_{i:08d} mock\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")

    print(f"Created {outfile} with {n_reads} reads")


if __name__ == "__main__":
    out_dir = Path(__file__).parent / "mock_fastq"
    out_dir.mkdir(exist_ok=True)

    make_mock_fastq(out_dir / "mock_euploid.fastq.gz", n_reads=2000, seed=1)
    make_mock_fastq(out_dir / "mock_t21.fastq.gz",     n_reads=2000, seed=2)
