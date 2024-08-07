from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="Click",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for CLICK analysis on PDB files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/Click",
    package_dir={"": "src"},
    packages=["Click"],  # Explicitly specify the package
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython",
        "pandas",
        "tqdm",
    ],
    entry_points={
        "console_scripts": [
            "click-analysis=Click.click:main",
        ],
    },
    package_data={
        "Click": ["bin/click"],
    },
)