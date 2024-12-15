# Matrix completion methods for pseudo-data generation 

This repo contains the code and data required to replicate the research of Matrix completion methods for pseudo-data generation of excess enthalpy data done by FS Middleton and JT Cripwell.
This work is in the process of being published for an article of the same name. 

## Table of Contents
1. [Installation and usage](#installation-and-usage)
3. [Folder Structure](#folder-structure)
4. [Contributing](#contributing)
5. [Contact](#contact)

## Installation and Usage

This repo can be cloned and used locally in a Julia editor for the UNIFAC predictions, and a 
Matlab editor for the Matlab code for array formations and experiments. [Julia](https://julialang.org/downloads/) can be installed locally, and [Maltab](https://matlab.mathworks.com/) can be accessed via a free online trial for 20 hours a month.
\\
This repo contains: 
* unifac:
    * The Julia code utilises the compounds expressed in the Excel spreadsheet to create UNIFAC (Do) predictions using the [Clapeyron library](https://github.com/ClapeyronThermo/Clapeyron.jl). 
    * This output is then loaded into Matlab. 
* data: All excess enthalpy data used in the research. The data is not present here due to permissions but an example is given in HEData.csv. Please contact Jamie Cripwell at cripwell@sun.ac.za for access to the data. 
* The src folder:
    * SVD_example: An example of singular value decomposition (SVD) of an array and how a scree plot is formed. 
    * ArrayFormation3way: The functions required to form a 3-way array from the excess enthalpy data provided, are used in ArrayFormation3way. 
    * Example_experiment: An example of an experiment conducted in the research.
    * The functions required to complete a matrix using the parallel completion method proposed in the paper and a simple SVD algorithm. This also includes the scripts for the initial filling guesses and calculation of the wSMSE, allowing the user to find the best rank for an array.

This repo should be used by:
* Creating UNIFAC predictions using UNIFACPredsJulia. Run the script in your IDE as it uses UNIFACParams.xlsx to run. Add the compounds you wish to evaluate to the sheet in the Excel sheet if you want to add any. This will output to `data/`.
* Transfer these predictions to a Matlab array using UNIFACClapeyron2Matlab in your Matlab editor.
* Create an array from experimental data stored in HEData.csv using ArrayFormation3way.
* Complete an example experiment using Example_experiment and the data you created. 

### Folder structure 

    . 
    ├── .github/                        # Github related files
    ├── data/                           # Original data as an Excel sheet. Used with ArrayFormation3way 
    ├── src/                            # Source code files
    │   ├── ArrayFormation3way.m        # Script for transforming the data into a usable sparse matrix
    │   ├── completion_2way_par.m       # Function for running one set of conditions for an experiment
    │   ├── Example_experiment.m        # An example of an experiment, using completion_2way_par
    │   ├── fill_data3.m                # Function filling an array with predictions of the entries
    │   ├── find_wmse_error.m           # Function to calculate the winsorized MSE error of predictions 
    │   ├── interp_data.m               # Function for interpolating data; aids in processing raw data 
    │   ├── missing_svd_par.m           # The matrix completion algorithm with thresholding
    │   └── SVD_example.m/              # Example of SVD used to complete a sparse matrix and draw scree plots
    ├── unifac/                         # Files relating to UNIFAC (Do) predictions using the [Clapeyron library](https://github.com/ClapeyronThermo/Clapeyron.jl)
    │   ├── UNIFACClapeyron2MATLAB.m    # Script for importing UNIFAC predictions into Matlab 
    │   ├── UNIFACParams.xlsx           # A file of the mixtures found in the experimental data and used in the Julia script
    │   └── UNIFACPredsJulia.jl         # Script creating UNIFAC (Do) predictions
    ├── citation.cff                    # How to cite this repo
    └── README.md                       # This document :) 



Please contact Francesca Middleton, @franmiddleton, or Jamie Cripwell, @jamiecripwell, for any queries.
