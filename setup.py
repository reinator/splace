from setuptools import setup, find_packages

setup(
    name="splace",
    version="4.0.0",
    description="SPLACE - Automation tool for viral and organellar phylogenomics",
    author="Renato Oliveira and Luan Rabelo",
    author_email="luan.rabelo@pq.itv.org",
    packages=find_packages(),
    scripts=['splace.py'],  # Make splace.py executable as a script
    install_requires=[
        "biopython",
        "requests",
        "pandas",
        "openpyxl",
        "syngenes",
        "pyyaml",
        "seaborn",
        "matplotlib",
        # 'mafft', 'trimal', 'iqtree' are binaries managed by conda, not pip
    ],
    python_requires=">=3.12",
)

