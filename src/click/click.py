"""
click.py

This module provides the Click class for performing CLICK analysis on PDB files
and generating pairwise alignments.

The Click class uses the CLICK tool to align protein structures and creates pairwise FASTA files.
"""

import os
import shutil
import subprocess
import tempfile
import logging
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from multiprocessing import Pool, cpu_count
import itertools
from typing import List, Tuple
from tqdm import tqdm
import pkg_resources

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

class ClickError(Exception):
    """Custom exception for CLICK related errors."""
    pass

class Click:
    """
    Performs CLICK analysis on PDB files and generates pairwise alignments.
    """

    def __init__(self, input_dir: Path, output_dir: Path, click_path: Path = None, cpu_percentage: float = 25):
        """
        Initialize the Click instance.

        Args:
            input_dir (Path): Directory containing input PDB files.
            output_dir (Path): Directory for output FASTA files.
            click_path (Path): Path to the CLICK executable. If None, defaults to the 'bin' directory.
            cpu_percentage (float): Percentage of CPU cores to use (default: 25).
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        if click_path is None:
            click_path = pkg_resources.resource_filename('click', 'bin/click')
        self.click_path = Path(click_path)
        self.cpu_cores = max(1, int(cpu_count() * (cpu_percentage / 100)))
        self.expected_fastas = 0
        self.actual_fastas = 0

    def run_analysis(self, pairwise_dir: Path = None) -> None:
        """
        Run the CLICK analysis pipeline to generate pairwise alignments.

        Args:
            pairwise_dir (Path): Directory to store pairwise FASTA files (default: None).
        """
        pairwise_dir = pairwise_dir or self.output_dir / "pairwise_fastas"
        pairwise_dir.mkdir(parents=True, exist_ok=True)

        pdb_files = list(self.input_dir.glob('*.pdb'))
        if not pdb_files:
            logging.warning(f"No PDB files found in {self.input_dir}")
            return

        self.expected_fastas = len(pdb_files) ** 2
        pairs = list(itertools.product(pdb_files, repeat=2))

        with Pool(processes=self.cpu_cores) as pool:
            list(tqdm(pool.imap_unordered(self._process_pair_wrapper,
                     [(pdb1, pdb2, pairwise_dir) for pdb1, pdb2 in pairs]),
                 total=len(pairs),
                 desc="Aligning pairs"))

        if self.expected_fastas != self.actual_fastas:
            logging.warning(f"Mismatch in FASTA file count. Expected: {self.expected_fastas}, Actual: {self.actual_fastas}")

    def _process_pair_wrapper(self, args):
        """
        Wrapper function for _process_pair to be used with Pool.imap_unordered.
        """
        return self._process_pair(*args)

    def _process_pair(self, pdb1: Path, pdb2: Path, pairwise_dir: Path) -> None:
        """
        Process a pair of PDB files with CLICK.

        Args:
            pdb1 (Path): Path to the first PDB file.
            pdb2 (Path): Path to the second PDB file.
            pairwise_dir (Path): Directory to store pairwise FASTA files.
        """
        try:
            self._run_click(pdb1, pdb2, pairwise_dir)
        except Exception as e:
            logging.error(f"Error processing {pdb1} and {pdb2}: {e}")

    def _run_click(self, pdb1: Path, pdb2: Path, pairwise_dir: Path) -> None:
        """
        Run CLICK on a pair of PDB files.

        Args:
            pdb1 (Path): Path to the first PDB file.
            pdb2 (Path): Path to the second PDB file.
            pairwise_dir (Path): Directory to store pairwise FASTA files.

        Raises:
            ClickError: If CLICK executable is not found or fails to run.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            pdb1_temp = os.path.join(temp_dir, os.path.basename(pdb1))
            pdb2_temp = os.path.join(temp_dir, os.path.basename(pdb2))
            shutil.copy(pdb1, pdb1_temp)
            shutil.copy(pdb2, pdb2_temp)

            pdb1_name = os.path.splitext(os.path.basename(pdb1))[0]
            pdb2_name = os.path.splitext(os.path.basename(pdb2))[0]
            alignment_output_prefix = f"{temp_dir}/{pdb1_name}-{pdb2_name}"

            command = [str(self.click_path), pdb1_temp, pdb2_temp]
            try:
                subprocess.run(command, check=True, capture_output=True, text=True, cwd=temp_dir)
            except subprocess.CalledProcessError as e:
                logging.error(f"CLICK failed for {pdb1} and {pdb2}: {e.stderr}")
                return
            except FileNotFoundError:
                raise ClickError(f"CLICK executable not found at {self.click_path}")

            expected_clique_file = f"{alignment_output_prefix}.pdb.1.clique"
            if not os.path.exists(expected_clique_file):
                raise ClickError(f"Clique file not found at expected location: {expected_clique_file}")

            df = self._parse_clique_file(expected_clique_file)
            seq_record1, seq_record2 = self._build_aligned_sequences(df, pdb1_name, pdb2_name)
            output_file = pairwise_dir / f"{pdb1_name}_{pdb2_name}.fasta"
            self._convert_to_fasta([seq_record1, seq_record2], output_file)
            self.actual_fastas += 1

    def _parse_clique_file(self, clique_file: str) -> pd.DataFrame:
        """Parse the CLICK clique file."""
        with open(clique_file, 'r') as file:
            header_line = next(i for i, line in enumerate(file) if line.strip().startswith('Chain'))

        df = pd.read_csv(clique_file, skiprows=header_line, sep=r'\s+',
                         names=['Chain1', 'ResNum1', 'ResName1', 'Atom1', 'Chain2', 'ResNum2', 'ResName2', 'Atom2'])

        return df[df['Atom1'] == 'CA'].astype({'ResNum1': int, 'ResNum2': int, 'ResName1': str, 'ResName2': str})

    def _build_aligned_sequences(self, df: pd.DataFrame, pdb1_name: str, pdb2_name: str) -> Tuple[SeqRecord, SeqRecord]:
        """Build aligned sequences from the parsed clique file."""
        seq1 = Seq(''.join(df['ResName1']))
        seq2 = Seq(''.join(df['ResName2']))
        return (SeqRecord(seq1, id=pdb1_name, description=""),
                SeqRecord(seq2, id=pdb2_name, description=""))

    def _convert_to_fasta(self, seq_records: List[SeqRecord], fasta_output_file: Path) -> None:
        """Convert sequence records to FASTA format and save to file."""
        with open(fasta_output_file, 'w') as outfile:
            for seq_record in seq_records:
                outfile.write(f">{seq_record.id}\n{''.join(seq_record.seq)}\n")

def click_analysis(input_dir: Path, output_dir: Path, click_path: Path = None, cpu_percentage: float = 25, pairwise_dir: Path = None) -> None:
    """
    Main function to run the CLICK analysis for generating pairwise alignments.

    Args:
        input_dir (Path): Directory containing input PDB files.
        output_dir (Path): Directory for output FASTA files.
        click_path (Path): Path to the CLICK executable (default: None).
        cpu_percentage (float): Percentage of CPU cores to use (default: 25).
        pairwise_dir (Path): Directory to store pairwise FASTA files (default: None).
    """
    analyzer = Click(input_dir, output_dir, click_path, cpu_percentage)
    analyzer.run_analysis(pairwise_dir)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Perform CLICK analysis and generate pairwise alignments")
    parser.add_argument("input_dir", type=Path, help="Directory containing input PDB files")
    parser.add_argument("output_dir", type=Path, help="Output directory for alignment results")
    parser.add_argument("--click_path", type=Path, help="Path to CLICK executable")
    parser.add_argument("--cpu_percentage", type=float, default=25, help="Percentage of CPU cores to use (default: 25)")
    parser.add_argument("--pairwise_dir", type=Path, help="Directory for pairwise alignments (if different from output_dir)")
    args = parser.parse_args()

    click_analysis(args.input_dir, args.output_dir, args.click_path, args.cpu_percentage, args.pairwise_dir)

if __name__ == "__main__":
    main()


