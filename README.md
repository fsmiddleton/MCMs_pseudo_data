# Matrix completion methods for pseudo-data generation 

This repository will contain all code, in MATLAB, used for the prediction of excess enthalpy data using matrix completion methods. 
This work is in the publication process for a an article of the same name 


## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
3. [Folder Structure](#folder-structure)
4. [Contributing](#contributing)
5. [Contact](#contact)

## Installation

Instructions on how to install and set up the project. Include dependencies, setup steps, and any configuration information.

## Usage

### Example:
git clone https://github.com/username/repository.git
cd repository
pip install -r requirements.txt

## Folder structure 
.
├── .github/                   # Raw and processed data
│   ├── raw/                # Original data, not modified
│   └── processed/          # Data that has been cleaned and prepared
├── src/                    # Source code files
│   ├── __init__.py         # Initialize the src module
│   ├── SVD_example.m       # Example of SVD used to complete a sparse matrix
│   └── utils/              # Utility functions
├── unifac/                  # Unit and integration tests
│   ├── UNIFACClapeyron2MATLAB.m    # Script for importing UNIFAC predictions into Matlab 
│   ├── UNIFACParams.xlsx           # A file containing the mixtures to predict excess enthalpy for using the Julia script
│   └── UNIFACPredsJulia.jl         # Script creating UNIFAC (Do) predictions
├── citation.cff            # How to cite this repo
├── requirements.txt        # Matlab and Julia dependencies
└── README.md               # Project readme

data/: Contains raw and processed data used in the project.
raw/: Original data, not modified.
processed/: Data that has been cleaned and prepared.
docs/: Contains documentation files, such as API references and guides.
src/: Contains the source code, including the main scripts and any utilities.
main.py: Main script to run the project.
utils/: Folder for utility functions used across the project.
tests/: Contains unit and integration tests for the project.
requirements.txt: Lists the Python dependencies required for the project.
README.md: The readme file you are currently reading.
