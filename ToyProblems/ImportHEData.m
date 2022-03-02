%% Importing excess enthalpy data into a tensor for use

%Francesca Middleton, 2022-03-02
clc
clear
%% Import the data of composition, component, temperature, real enthalpy and excess enthalpy
data = readtable('HEData.xlsx','ReadVariableNames',true); % change sheet to include certain functional groups as the main site of data collection 

comp = table2array(data(:,6));
temp = table2array(data(:,7));
HE  = table2array(data(:,12));

%% Interpolate data 
% at one temp (ignore pressure) for each mixture to find values at certain concentrations 
% Find all data points with the same temp, FuncGroup1, FuncGroup2,
% ChainLength1 and ChainLength2 and populate comp_temp with compositions
% and HE_temp with HE data 

%collect first data point 
%find all data points that match this
interp_index = zeros(length(data.FunctionalGroup1),1);%variable to save the indexes that have been interpolated
index = 1;%current index, 
%collect first data point
% temp.functional_group.one = data.FunctionalGroup1(index);
% temp.functional_group.two = data.FunctionalGroup2(index);
% temp.chain_length.one = data.ChainLength1(index);
% temp.chain_length.two = data.ChainLength2(index);

% Specify the mixtures wanted in the matrix. The algorithm will find all
% combinations of functional group 1 and 2. Only organic molecules were considered here  
func_groups.one = ['Alkane', 'Primaryalcohol'];% 'Secondaryalcohol', 'Ketone', 'Alkene','Cycloalkane'];
func_groups.two = ['Alkane', 'Primaryalcohol'];
max_chain_length = 10; 
T = 298.15; % temperature in kelvin used for the matrix
% P = 101.33; % pressure in kPa used for the matrix

conc_interval = 0:0.1:1;
%First restrict data to relevant temperatures. Can add a for loop here
%later 
T_loc = find(data.Temperature==T);
temp_data = data(T_loc, :);

for func1= func_groups.one
    % restrict to mixtures containing the wanted functional group in
    % position 1
    func1_loc = find(temp_data.FunctionalGroup1==func1); % you can't directly compare strings 
    temp1 = temp_data(func1_loc);
    for func2= func_groups.two
        % restrict to mixtures containing the wanted functional group in
        % position 2
        func2_loc = find(temp1.FunctionalGroup2==func2);
        temp2 = temp1(func2_loc);
        % now loop through chain lengths 
        for i =1:max_chain_length
            chain1_loc = find(temp2.ChainLength1==i);
            temp3 = temp2(chain1_loc)
            for j =1:max_chain_length
                chain2_loc = find(temp3.ChainLength2==j);
                temp4 = temp3(chain2_loc); % this is the data we interpolate 
                [HE(i), uncertainty(i)]=interp_data(temp4, conc_interval);
            end 
        end 
        
    end 
end 


function [HE, uncertainty]=interp_data(data, conc_interval)
    % Interpolate data for use in the matrices
    % Inputs specify the mixture to be captured and the concentration
    % interval for which to create interpolate points 
    % data = the whole table of 
    % func_group = functional groups 1 and 2 in func_group.1 and ...2
    % chain_length = chain lengths of groups 1 and 2 in chain_length.1 and .2 
    % T = temperature in Kelvin 
    % conc_interval = list of concentrations for interpolation 
    
    % Ouputs is the data produced by the interpolation, ready for
    % matrix populating 
    % HE = excess enthalpy data at each concentration value for the mixture
    % uncertainty = uncertainty associated with each data point 
    
    HE = zeros(length(conc),1);
    for i=1:length(conc)
        %populate the data 
        HE(i) = interp_data(conc(i));
    end 

end 
