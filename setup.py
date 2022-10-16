from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
	name="girafr",  
	version="1.0.7", 
	description="A package to detect mutations and CRISPR-Cas9 editing effects in single cell CRISPR screens", 
	long_description=long_description, 
	long_description_content_type="text/markdown",  
	url="https://github.com/FunctionalEpigeneticsLab/GiRAFR", 
	author="Laboratory for Functional Epigentics, KU Leuven",  
	author_email="qian.yu@kuleuven.be", 
	classifiers=[  
		"Development Status :: 3 - Alpha",
		"Topic :: Scientific/Engineering :: Bio-Informatics",
		"License :: OSI Approved :: MIT License",
		"Programming Language :: Python :: 3",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
		"Programming Language :: Python :: 3.10",
		"Programming Language :: Python :: 3 :: Only",
	],
	keywords="gRNA, mutation, single cell CRISPR screen", 
	packages=find_packages(),
	python_requires=">=3.7, <4",
	install_requires=[
		"pysam",
		"Bio",
		"sklearn",
		"pandas >=1.0.0, ==1.*",
		"numpy"],  
	entry_points={
		"console_scripts": [
			"girafr = girafr.__main__:main",
		]
	},
	project_urls={  
		"Bug Reports": "https://github.com/FunctionalEpigeneticsLab/GiRAFR/issues",
		"Source": "https://github.com/FunctionalEpigeneticsLab/GiRAFR",
	},
)
