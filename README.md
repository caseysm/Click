# CLICK Analysis Package

This package provides tools for performing CLICK analysis on PDB files and generating pairwise alignments.

## Installation

To install this package directly from GitHub, use the following pip command:

```
pip install git+https://github.com/yourusername/click.git
```

Replace `yourusername` with your actual GitHub username.

## Usage

After installation, you can use the package in your Python scripts:

```python
from click import ClickAnalyzer

analyzer = ClickAnalyzer(input_dir="path/to/pdb/files", output_dir="path/to/output")
analyzer.run_analysis()
```

You can also use the command-line interface:

```
click-analysis path/to/pdb/files path/to/output
```

For more options, use:

```
click-analysis --help
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.