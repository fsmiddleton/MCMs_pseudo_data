%% Importing UNIFAC predictions 
% Imports UNIFAC predictions made in Clapeyron to MATLAB and writes into a
% MATLAB file 

%% Import data 
Temps = [243.15; 253.15; 263.15; 273.15; 283.15; 288.15; 290.15; 293.15; 296.15; 298.15; 303.15;307.5;
    309.5; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15; 363.15];
func_groups= {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 
    'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
%import compound data 
compounds = readtable("UNIFACParams.xlsx", 'Sheet', "compounds");
functionalgroup = table2cell(compounds(:,1));
chainlength = table2array(compounds(:,2));
compounds = table2array(compounds(:,3));

conc_interval = 0:0.01:1;
he = zeros(length(conc_interval),size(data,1));

for T = Temps' % for each excel file 
    %import data 
    data = readtable(strcat("UNIFACPreds", num2str(T), ".xlsx"));
    mixturesIn = table2array(data(:,1:2));% cell array
    %declare important variables 
    mixture = zeros(4,size(data,1));
    he = zeros(length(conc_interval),size(data,1));
    for i = 1:size(data,1) % for each mixture - populate the mixture and excess enthalpy arrays 
        % find the mixture
        comp1 = mixturesIn(i,1);
        comp2 = mixturesIn(i,2);
        ind1 = find(strcmp(comp1,compounds));
        ind2 = find(strcmp(comp2,compounds));
        %functional groups - the index is the number of the functional
        %group
        mixture(1,i) = find(strcmp(functionalgroup(ind1), func_groups));
        mixture(3,i) = find(strcmp(functionalgroup(ind2), func_groups));
        %chain lengths
        mixture(2,i) = chainlength(ind1);
        mixture(4,i) = chainlength(ind2);
        %excess enthalpy
        he(:,i) = table2array(data(i,3:end))';
    end 
    %export new file 
    filename = strcat('heUNIFACforT=', num2str(T),'.mat');
    save(filename)
end 
%%
clc
clear
T = 243.15;
conc_interval = 0:0.01:1;
data = readtable(strcat("UNIFACPreds", num2str(T), ".xlsx"));
%data = table2array(data);
mixturesIn = data(:,1:2);
mixture = zeros(4,size(data,1));
%%
data = table2array(data);
mixture = zeros(4,size(data,1));
conc_interval = 0:0.01:1;
func_groups= {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
func_groups.two = {'Alkane', 'Primaryalcohol'};%, 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
compounds = readtable("UNIFACParams.xlsx", 'Sheet', "compounds");
%%
i=1;
comp1 = mixturesIn(i,1);
comp2 = mixturesIn(i,2);
compounds = readtable("UNIFACParams.xlsx", 'Sheet', "compounds");
functionalgroup = table2cell(compounds(:,1));
chainlength = table2array(compounds(:,2));
compounds = table2array(compounds(:,3));
ind1 = find(strcmp(comp1,compounds(:,1)));
ind2 = find(strcmp(comp2,compounds(:,1)));
%functional groups
mixture(1,1) = find(strcmp(table2cell(compoundkey(ind1,1)), func_groups));
mixture(3,1) = find(strcmp(table2cell(compoundkey(ind2,1)), func_groups));
%chain lengths
mixture(2,1) = chainlength(ind1);
mixture(4,1) = chainlength(ind2);