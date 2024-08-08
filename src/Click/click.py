import os
import shutil
import subprocess
import tempfile
import logging
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from multiprocessing import Pool, cpu_count
import itertools
from typing import List, Tuple, Optional
from tqdm import tqdm
import pkg_resources

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

class ClickError(Exception):
    """Custom exception for CLICK related errors."""
    pass

class Click:
    """
    Performs CLICK analysis on PDB files and generates pairwise alignments and a similarity matrix.
    """

    def __init__(self, input_dir: Path, output_dir: Path, click_path: Path = None):
        """
        Initialize the Click instance.

        Args:
            input_dir (Path): Directory containing input PDB files.
            output_dir (Path): Directory for output files.
            click_path (Path): Path to the CLICK executable. If None, defaults to the 'bin' directory.
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        if click_path is None:
            click_path = pkg_resources.resource_filename('Click', 'bin/click')
        self.click_path = Path(click_path)
        self.expected_fastas = 0
        self.actual_fastas = 0
        self.structure_overlaps = {}
        self.pdb_names = []

    def pdb_vs_pdb(self, pdb1: Path, pdb2: Path, output_overlap: bool = True, output_clique: bool = False, output_aligned_pdb: bool = False) -> float:
        """
        Compare two PDB files using CLICK.

        Args:
            pdb1 (Path): Path to the first PDB file.
            pdb2 (Path): Path to the second PDB file.
            output_overlap (bool): Whether to output the overlap score (default: True).
            output_clique (bool): Whether to output the clique file (default: False).
            output_aligned_pdb (bool): Whether to output the aligned PDB files (default: False).

        Returns:
            float: The structure overlap score.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            pdb1_temp = Path(temp_dir) / pdb1.name
            pdb2_temp = Path(temp_dir) / pdb2.name
            shutil.copy(pdb1, pdb1_temp)
            shutil.copy(pdb2, pdb2_temp)

            pdb1_name = pdb1.stem
            pdb2_name = pdb2.stem
            alignment_output_prefix = f"{temp_dir}/{pdb1_name}-{pdb2_name}"

            command = [str(self.click_path), str(pdb1_temp), str(pdb2_temp)]
            try:
                subprocess.run(command, check=True, capture_output=True, text=True, cwd=temp_dir)
            except subprocess.CalledProcessError as e:
                logging.error(f"CLICK failed for {pdb1} and {pdb2}: {e.stderr}")
                return 0.0
            except FileNotFoundError:
                raise ClickError(f"CLICK executable not found at {self.click_path}")

            clique_file = Path(f"{alignment_output_prefix}.pdb.1.clique")
            if not clique_file.exists():
                raise ClickError(f"Clique file not found at expected location: {clique_file}")

            overlap = self._parse_overlap_from_clique(clique_file)

            if output_overlap:
                print(f"Structure overlap between {pdb1_name} and {pdb2_name}: {overlap:.2f}")

            if output_clique:
                shutil.copy(clique_file, self.output_dir / clique_file.name)

            if output_aligned_pdb:
                aligned_pdb_files = list(Path(temp_dir).glob(f"{pdb1_name}-{pdb2_name}.pdb.*"))
                for aligned_pdb in aligned_pdb_files:
                    shutil.copy(aligned_pdb, self.output_dir / aligned_pdb.name)

            return overlap

    def all_vs_all(self, cpu_percentage: float = 25, output_similarity_matrix: bool = True, output_clique: bool = False, output_aligned_pdb: bool = False) -> Optional[np.ndarray]:
        """
        Perform all-vs-all comparison of PDB files in the input directory.

        Args:
            cpu_percentage (float): Percentage of CPU cores to use (default: 25).
            output_similarity_matrix (bool): Whether to output the similarity matrix (default: True).
            output_clique (bool): Whether to output the clique files (default: False).
            output_aligned_pdb (bool): Whether to output the aligned PDB files (default: False).

        Returns:
            Optional[np.ndarray]: The similarity matrix if output_similarity_matrix is True, otherwise None.
        """
        pdb_files = list(self.input_dir.glob('*.pdb'))
        if not pdb_files:
            logging.warning(f"No PDB files found in {self.input_dir}")
            return None

        self.pdb_names = [pdb.stem for pdb in pdb_files]
        self.expected_fastas = len(pdb_files) ** 2
        pairs = list(itertools.product(pdb_files, repeat=2))

        cpu_cores = max(1, int(cpu_count() * (cpu_percentage / 100)))
        with Pool(processes=cpu_cores) as pool:
            list(tqdm(pool.imap_unordered(self._process_pair_wrapper,
                     [(pdb1, pdb2, output_clique, output_aligned_pdb) for pdb1, pdb2 in pairs]),
                 total=len(pairs),
                 desc="Aligning pairs"))

        if self.expected_fastas != self.actual_fastas:
            logging.warning(f"Mismatch in FASTA file count. Expected: {self.expected_fastas}, Actual: {self.actual_fastas}")

        sim_matrix = self._calculate_similarity_matrix()

        if output_similarity_matrix:
            self._save_similarity_matrix(sim_matrix)
            return sim_matrix
        else:
            return None

    def _process_pair_wrapper(self, args):
        """Wrapper function for _process_pair to be used with Pool.imap_unordered."""
        return self._process_pair(*args)

    def _process_pair(self, pdb1: Path, pdb2: Path, output_clique: bool, output_aligned_pdb: bool) -> None:
        """Process a pair of PDB files with CLICK."""
        try:
            self._run_click(pdb1, pdb2, output_clique, output_aligned_pdb)
        except Exception as e:
            logging.error(f"Error processing {pdb1} and {pdb2}: {e}")

    def _run_click(self, pdb1: Path, pdb2: Path, output_clique: bool, output_aligned_pdb: bool) -> None:
        """Run CLICK on a pair of PDB files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            pdb1_temp = Path(temp_dir) / pdb1.name
            pdb2_temp = Path(temp_dir) / pdb2.name
            shutil.copy(pdb1, pdb1_temp)
            shutil.copy(pdb2, pdb2_temp)

            pdb1_name = pdb1.stem
            pdb2_name = pdb2.stem
            alignment_output_prefix = f"{temp_dir}/{pdb1_name}-{pdb2_name}"

            command = [str(self.click_path), str(pdb1_temp), str(pdb2_temp)]
            try:
                subprocess.run(command, check=True, capture_output=True, text=True, cwd=temp_dir)
            except subprocess.CalledProcessError as e:
                logging.error(f"CLICK failed for {pdb1} and {pdb2}: {e.stderr}")
                return
            except FileNotFoundError:
                raise ClickError(f"CLICK executable not found at {self.click_path}")

            clique_file = Path(f"{alignment_output_prefix}.pdb.1.clique")
            if not clique_file.exists():
                raise ClickError(f"Clique file not found at expected location: {clique_file}")

            self._parse_overlap_from_clique(clique_file, pdb1_name, pdb2_name)

            if output_clique:
                shutil.copy(clique_file, self.output_dir / clique_file.name)

            if output_aligned_pdb:
                aligned_pdb_files = list(Path(temp_dir).glob(f"{pdb1_name}-{pdb2_name}.pdb.*"))
                for aligned_pdb in aligned_pdb_files:
                    shutil.copy(aligned_pdb, self.output_dir / aligned_pdb.name)

            self.actual_fastas += 1

    def _parse_overlap_from_clique(self, clique_file: Path, pdb1_name: Optional[str] = None, pdb2_name: Optional[str] = None) -> float:
        """Parse the structure overlap from the clique file."""
        with open(clique_file, 'r') as file:
            for line in file:
                if line.startswith("Structure Overlap"):
                    overlap = float(line.split('=')[1].strip()) / 100.0
                    if pdb1_name and pdb2_name:
                        self.structure_overlaps[(pdb1_name, pdb2_name)] = overlap
                        self.structure_overlaps[(pdb2_name, pdb1_name)] = overlap
                    return overlap
        raise ClickError(f"Structure Overlap not found in clique file: {clique_file}")

    def _calculate_similarity_matrix(self) -> np.ndarray:
        """Calculate similarity matrix from structure overlap scores."""
        n = len(self.pdb_names)
        sim_matrix = np.zeros((n, n))

        for i, pdb1 in enumerate(self.pdb_names):
            for j, pdb2 in enumerate(self.pdb_names):
                if i == j:
                    sim_matrix[i, j] = 1.0  # Self-similarity
                else:
                    overlap = self.structure_overlaps.get((pdb1, pdb2), 0)
                    sim_matrix[i, j] = overlap

        return sim_matrix

    def _save_similarity_matrix(self, sim_matrix: np.ndarray) -> None:
        """Save the similarity matrix to a CSV file."""
        df = pd.DataFrame(sim_matrix, index=self.pdb_names, columns=self.pdb_names)
        output_file = self.output_dir / "similarity_matrix.csv"
        df.to_csv(output_file)
        logging.info(f"Similarity matrix saved to {output_file}")

def click_analysis(input_dir: Path, output_dir: Path, click_path: Path = None, cpu_percentage: float = 25,
                   mode: str = "all_vs_all", pdb1: Path = None, pdb2: Path = None,
                   output_similarity_matrix: bool = True, output_overlap: bool = True,
                   output_clique: bool = False, output_aligned_pdb: bool = False) -> None:
    """
    Main function to run the CLICK analysis for generating pairwise alignments and similarity matrix.

    Args:
        input_dir (Path): Directory containing input PDB files.
        output_dir (Path): Directory for output files.
        click_path (Path): Path to the CLICK executable (default: None).
        cpu_percentage (float): Percentage of CPU cores to use for all_vs_all mode (default: 25).
        mode (str): Analysis mode, either "all_vs_all" or "pdb_vs_pdb" (default: "all_vs_all").
        pdb1 (Path): Path to the first PDB file for pdb_vs_pdb mode (default: None).
        pdb2 (Path): Path to the second PDB file for pdb_vs_pdb mode (default: None).
        output_similarity_matrix (bool): Whether to output the similarity matrix in all_vs_all mode (default: True).
        output_overlap (bool): Whether to output the overlap score in pdb_vs_pdb mode (default: True).
        output_clique (bool): Whether to output the clique file(s) (default: False).
        output_aligned_pdb (bool): Whether to output the aligned PDB files (default: False).
    """
    analyzer = Click(input_dir, output_dir, click_path)

    if mode == "all_vs_all":
        analyzer.all_vs_all(cpu_percentage, output_similarity_matrix, output_clique, output_aligned_pdb)
    elif mode == "pdb_vs_pdb":
        if pdb1 is None or pdb2 is None:
            raise ValueError("Both pdb1 and pdb2 must be provided for pdb_vs_pdb mode")
        analyzer.pdb_vs_pdb(pdb1, pdb2, output_overlap, output_clique, output_aligned_pdb)
    else:
        raise ValueError("Invalid mode. Choose either 'all_vs_all' or 'pdb_vs_pdb'")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Perform CLICK analysis on PDB files")

    # Common arguments
    parser.add_argument("--click_path", type=Path, help="Path to CLICK executable")
    parser.add_argument("--output_clique", action="store_true", help="Output clique file(s)")
    parser.add_argument("--output_aligned_pdb", action="store_true", help="Output aligned PDB files")
    parser.add_argument("--mode", choices=["all_vs_all", "pdb_vs_pdb"], required=True, help="Analysis mode")

    # All-vs-all mode arguments
    parser.add_argument("--input_dir", type=Path, help="Directory containing input PDB files (for all_vs_all mode)")
    parser.add_argument("--output_dir", type=Path, help="Output directory for results (for all_vs_all mode)")
    parser.add_argument("--cpu_percentage", type=float, default=25,
                        help="Percentage of CPU cores to use for all_vs_all mode (default: 25)")
    parser.add_argument("--no_similarity_matrix", action="store_true",
                        help="Don't output similarity matrix in all_vs_all mode")

    # PDB-vs-PDB mode arguments
    parser.add_argument("--pdb1", type=Path, help="Path to first PDB file (for pdb_vs_pdb mode)")
    parser.add_argument("--pdb2", type=Path, help="Path to second PDB file (for pdb_vs_pdb mode)")
    parser.add_argument("--output_fasta", type=Path, help="Output FASTA file path (for pdb_vs_pdb mode)")
    parser.add_argument("--no_overlap", action="store_true", help="Don't output overlap score in pdb_vs_pdb mode")

    args = parser.parse_args()

    if args.mode == "all_vs_all":
        if not args.input_dir or not args.output_dir:
            parser.error("--input_dir and --output_dir are required for all_vs_all mode")
        click_analysis(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            click_path=args.click_path,
            cpu_percentage=args.cpu_percentage,
            mode="all_vs_all",
            output_similarity_matrix=not args.no_similarity_matrix,
            output_clique=args.output_clique,
            output_aligned_pdb=args.output_aligned_pdb
        )
    elif args.mode == "pdb_vs_pdb":
        if not args.pdb1 or not args.pdb2 or not args.output_fasta:
            parser.error("--pdb1, --pdb2, and --output_fasta are required for pdb_vs_pdb mode")
        click_analysis(
            input_dir=Path(args.pdb1).parent,  # Use parent directory of pdb1 as input_dir
            output_dir=Path(args.output_fasta).parent,  # Use parent directory of output_fasta as output_dir
            click_path=args.click_path,
            mode="pdb_vs_pdb",
            pdb1=args.pdb1,
            pdb2=args.pdb2,
            output_fasta=args.output_fasta,
            output_overlap=not args.no_overlap,
            output_clique=args.output_clique,
            output_aligned_pdb=args.output_aligned_pdb
        )


if __name__ == "__main__":
    main()