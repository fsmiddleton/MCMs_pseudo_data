# Matrix completion methods for pseudo-data generation 

This repo contains the code and data required to replicate the research of Matrix completion methods for pseudo-data generation of excess enthalpy data done by FS Middleton and JT Cripwell.
This work is in the process of being published for an article of the same name. 

## Table of Contents
1. [Installation and usage](#installation-and-usage)
2. [Folder Structure](#folder-structure)
3. [Contact](#contact)

## Installation and Usage

This repo can be cloned and used locally in a Julia editor for the UNIFAC predictions, and a 
Matlab editor for the Matlab code for array formations and experiments. [Julia](https://julialang.org/downloads/) can be installed locally, and [Maltab](https://matlab.mathworks.com/) can be accessed via a free online trial for 20 hours a month.
\
This repo should be used by:
* Creating UNIFAC predictions using `UNIFACPredsJulia.jl`. Run the script in your IDE as it uses `UNIFACParams.xlsx` to run. Add the compounds you wish to evaluate to the sheet in the Excel sheet if you want to add any. This will output to `data/`.
* Transfer these predictions to a Matlab array using `UNIFACClapeyron2Matlab.m` in your Matlab editor. Ensure the temperatures you use and the functional groups are the same e are consistent with the contents of `UNIFACParams.xlsx`.
* Create an array from experimental data stored in HEData.csv using `ArrayFormation3way.m`. 
* Complete an example experiment using `Example_experiment.m` and the data you created. Ensure the filename you outputted in `ArrayFormation3way.m` is the same filename as the filename you import in this script. Adjust the parameters as desired: r, T, fillmethod, maxiter, and thresholdperc.   
\
This repo contains: 
* unifac:
    * `UNIFACPredsJulia.jl` utilises the compounds expressed in the Excel spreadsheet to create UNIFAC (Do) predictions using the [Clapeyron library](https://github.com/ClapeyronThermo/Clapeyron.jl). 
    * This output is then loaded into Matlab using `UNIFACClapeyron2MATLAB.m`. 
* data: Some mock excess enthalpy data. The data used in the research is not present here due to permissions but an example is given in `HEData.csv`. Please contact Jamie Cripwell at cripwell@sun.ac.za for access to the data. 
* The src folder:
    * `SVD_example.m`: An example of singular value decomposition (SVD) of an array and how a scree plot is formed. 
    * `ArrayFormation3way.m`: A 3-way array is formed from data using this script. This calls `interp_data.m`.
    * `Example_experiment.m`: An example of an experiment conducted in the research. `completion_2way_par.m` is called here to perform array completion, which in turn calls `missing_svd_par.m`. `find_wmse_error.m` is used to find error metrics. 

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


### Contact 
Please contact Francesca Middleton, @franmiddleton, or Jamie Cripwell, @jamiecripwell, for any queries.
