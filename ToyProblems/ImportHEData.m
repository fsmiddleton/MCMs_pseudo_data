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
func_groups.one = {'Alkane', 'Primaryalcohol'};% 'Secondaryalcohol', 'Ketone', 'Alkene','Cycloalkane'];
func_groups.two = { 'Primaryalcohol'};
max_chain_length = 10; 
T = 298.15; % temperature in kelvin used for the matrix
% a moderate allowance for different experimental values 
Tupper = T+1;
Tlower = T-1;
% P = 101.33; % pressure in kPa used for the matrix
% pressure is ignored due to very small variation with pressure of HE

conc_interval = 0.1:0.1:0.9;
c = length(conc_interval);
%First restrict data to relevant temperatures. Can add a for loop here
%later 

T_loc = find(data.Temperature<Tupper & data.Temperature>Tlower);
temp_data = data(T_loc, :);

% create matrices to populate with each mixture's data 
HE_data = nan(length(conc_interval), (length(func_groups.one)*length(func_groups.two)*max_chain_length^2));
dim1= 40;
conc_original = zeros(dim1, (length(func_groups.one)*length(func_groups.two)*max_chain_length^2));
HE_original = zeros(dim1, (length(func_groups.one)*length(func_groups.two)*max_chain_length^2));
mixture = zeros(4, (length(func_groups.one)*length(func_groups.two)*max_chain_length^2));
f1=0;
ind=0;
for func1= func_groups.one
    f1=f1+1;
    f2=0;
    % restrict to mixtures containing the wanted functional group in
    % position 1
    func1_loc= find(strcmp(temp_data.FunctionalGroup1, func1));
    temp1 = temp_data(func1_loc,:);
    for func2= func_groups.two
        f2=f2+1;
        % restrict to mixtures containing the wanted functional group in
        % position 2
        func2_loc = find(strcmp(temp1.FunctionalGroup2, func2));
        temp2 = temp1(func2_loc,:);
        % now loop through chain lengths 
        for i =1:max_chain_length
            chain1_loc = find(temp2.Chainlength1==i);
            temp3 = temp2(chain1_loc,:);
            for j =1:max_chain_length
                ind =ind+1;% this index is wrong
                mixture(1,ind)= f1;
                mixture(2,ind)= f2;
                mixture(3,ind)= i;
                mixture(4,ind)= j;
                if i==j && f1==f2
                    %excess enthalpy of a mixture of the same two pure
                    %components = 0
                    HE_data(:,ind) = zeros(length(conc_interval),1);
                else 
                    chain2_loc = find(temp3.Chainlength2==j);
                    temp4 = temp3(chain2_loc,:); % this is the data we interpolate 
                    if size(temp4,1)>1
                        % interpolate data and populate matrix of all data
                        [HE_data(:,ind), uncertainty(:,f2*f1*i*j)]=interp_data(temp4, conc_interval);
                        % save original data in the same order
                        conc_original(1:size(temp4,1),ind) = temp4.Compositioncomponent1;
                        HE_original(1:size(temp4,1),ind) = temp4.Excessenthalpy;
                    end 
                end 
            end 
            
        end 
        
    end 
end
%% Check the interpolation worked nicely 
clf
mix =13;
plot(conc_original(:,mix), HE_original(:,mix), 'bo')
hold on 
plot(conc_interval, HE_data(:,mix), 'LineWidth', 1) 

legend('Experimental', 'Interpolated', 'Location', 'northwest')

xlabel('Composition component 1 (mol/mol)')
ylabel('Excess enthalpy (kJ/mol)')
disp(mixture(:,mix))
%% Populate matrices
dim1 = max_chain_length*length(func_groups.one);
dim2 = max_chain_length*length(func_groups.two);
dim3 = length(conc_interval);
HE_matrix= zeros(dim1,dim2,dim3);
row = 0;
for c = 1:length(conc_interval) % for each concentration in the interval 
    % not sure if HE_data reshaped nicely 
    for j = 0:(length(func_groups.one)-1)
        HE_matrix(1+dim1/length(func_groups.one)*j:dim1/length(func_groups.one)*(j+1),:,c) = reshape(HE_data(c,1+max_chain_length*dim2*j:max_chain_length*dim2*(j+1)),[dim1/length(func_groups.one),dim2]);
    end 
end 

%find missing rows and columns 
missing.ind = find(isnan(HE_matrix));
% check for nan rows and column 
[missing.i, missing.j] = find(isnan(HE_matrix));
% remove all data that is nan
% remove from mixtures, HE_matrix, HE_original, comp_original 

% rows 
B=HE_matrix(sum(isnan(HE_matrix),2)==0);
%% Export to excel spreadsheet 
filename = 'HEMatrix298.15.2.xlsx';
%create table with all this information
for i = 1:length(conc_interval)
    T = array2table(HE_matrix(:,:,i));
    writetable(T,filename,'Sheet',num2str(conc_interval(i)))
end 

%%
function [HE, uncertainty]=interp_data(data, conc_interval)
    % Interpolate data for use in the matrices
    % Inputs specify the mixture to be captured and the concentration
    % interval for which to create interpolate points
    
    % data = the whole table of data to be interpolated
    % conc_interval = list of concentrations for interpolation 
    
    % Ouputs is the data produced by the interpolation, ready for
    % matrix populating 
    % HE = excess enthalpy data at each concentration value for the mixture
    % uncertainty = uncertainty associated with each data point, to be
    % added
    HE_original = data.Excessenthalpy;
    comp = data.Compositioncomponent1;
    % remove the 0 and 1 data points 
    ind_keep = find(HE_original ~=0);
    
    HE_new= HE_original(ind_keep);
    comp_new = comp(ind_keep);
    %interpolate the data 
    
    [p,S] = polyfit(comp_new, HE_new,7);
    [HE,uncertainty] = polyval(p,conc_interval, S);
    %populate the data 

end 
