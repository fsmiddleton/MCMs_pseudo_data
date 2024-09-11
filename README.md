# Matrix completion methods for pseudo-data generation 

This repository will contain all code, in MATLAB, used for the prediction of excess enthalpy data using matrix completion methods. 
This work is in the publication process for a an article of the same name 


## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
3. [Folder Structure](#folder-structure)
4. [Contributing](#contributing)
5. [Contact](#contact)

## Installation and Usage

This repo can be cloned and used locally in a Julia editor for the UNIFAC predictions, and a 
Matlab editor for the Matlab code for array formations and experiments. 

## Usage

This repo contains the code and data required to replicate the article titled Matrix completion methods for pseudo-data generation 
from FS Middleton and JT Cripwell.
\\
This repo contains: 
* unifac:
    * The Julia code which utilises the compounds expressed in the Excel spreadsheet to create UNIFAC (Do) predictions using the Clapeyron library. 
    * This output is then loaded into Matlab. 
* The src folder:
    * An example of singular value decomposition (SVD) of an array and how a scree plot is formed. 
    * The functions required to perform the formation of a 3-way array from the excess enthalpy data provided. 
    * The functions required to complete a matrix using the parallel completion method proposed in the paper and a simple SVD algorithm. This also includes the scripts for the initial filling guesses and calculation of the wSMSE, allowing the user to find the best rank for an array.

## Folder structure 
.
├── .github/                        # Github related files
├── data/                           # Original data, not modified
├── src/                            # Source code files
│   ├── ArrayFormation3way.m        # How the data in data/ can be transformed into a usable sparse matrix
│   ├── completion_2way_par.m       # The completion algorithm 
│   ├── Example_experiment.m        # An example of an experiment, which can be tailored to any of the data
│   ├── fill_data3.m                # Fills an array with predictions of the entries
│   ├── find_wmse_error.m           # Function to calculate the winsorized MSE error of predictions 
│   ├── interp_data.m               # Interpolates data and aids in processing raw data 
│   ├── missing_svd_par.m           # The completion algorithm 
│   └── SVD_example.m/              # Example of SVD used to complete a sparse matrix and draw scree plots
├── unifac/                         # Files relating to UNIFAC (Do) predictions using  the [Clapeyron library](https://github.com/ClapeyronThermo/Clapeyron.jl)
│   ├── UNIFACClapeyron2MATLAB.m    # Script for importing UNIFAC predictions into Matlab 
│   ├── UNIFACParams.xlsx           # A file containing the mixtures to predict excess enthalpy for using the Julia script
│   └── UNIFACPredsJulia.jl         # Script creating UNIFAC (Do) predictions
├── citation.cff                    # How to cite this repo
└── README.md                       # Project readme

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
