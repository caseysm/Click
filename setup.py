from setuptools import setup, find_packages
from setuptools.command.install import install
import os
import stat

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        # Set executable permissions for the CLICK binary
        click_path = os.path.join(self.install_lib, 'Click', 'bin', 'click')
        if os.path.exists(click_path):
            st = os.stat(click_path)
            os.chmod(click_path, st.st_mode | stat.S_IEXEC)

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
    packages=find_packages(where="src"),
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
        "requests",
    ],
    entry_points={
        "console_scripts": [
            "click=Click.click:main",
        ],
    },
    cmdclass={
        'install': PostInstallCommand,
    },
    package_data={
        "Click": ["bin/click"],
    },
)