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
%% Find errors of all UNIFAC (Do) predictions 
clc
clear
%initiate necessary vars 
Temps = [ 283.15, 288.15, 293.15, 298.15, 303.15, 307.5,309.5, 313.15, 318.15, 323.15, 363.15];
smse_overall = zeros(length(Temps),1);
nopreds = zeros(length(Temps),1);
AARD_overall = zeros(length(Temps),1);
count=0;% temperature counter 
for T = Temps% for each excel file
    count=count+1;
    %import experimental data 
    filename = strcat('HEData3wayPolyAll',num2str(T),'.mat');
    load(filename)
    %columns = mixtures throughout
    mixture_exp = mixture;
    he_exp = HE_data; 
    conc_exp = conc_interval;
    %import data for unifac predictions
    filename = strcat('heUNIFACforT=', num2str(T),'.mat');
    load(filename)
    % save predictions 
    mixture_pred = mixture;
    he_pred = zeros(size(HE_data)); % same size as the experimental data 
    conc_pred = conc_interval;
    conc_indices = 11:10:91; %known from previous code - indices of predictions 
    ind_keep = 1:size(mixture_exp,2);
    for i = 1:size(mixture_exp,2)
        [~,mix_ind] = (ismember(mixture_exp(:,i)',mixture_pred','rows'));
        
        if mix_ind ==0
            tempmix = [mixture_exp(3:4,i); mixture_exp(1:2,i)];
            [~,mix_ind] = (ismember(tempmix', mixture_pred', 'rows'));
        end 
        if mix_ind ==0 
            disp('no prediction')
            disp(mixture_exp(:,i))
            ind_keep(i) = 0;
        else
            he_pred(:,i) = he(conc_indices,mix_ind);
        end
        if he_pred(5,i) ==0 % prediction will ruin overall errors
            ind_keep(i) =0;
        end 
    end 
    
    he_exp = he_exp(:,find(ind_keep));
    he_pred = he_pred(:,find(ind_keep));
    
    err = he_exp-he_pred;
    
    %remove zeros 
    smse_sys = sqrt(sum(err.^2)/size(he_pred,1)); % per system
    smse_overall(count) = sqrt(sum(sum(err.^2))/prod(size(he_pred)));
    nopreds(count) = size(he_pred,1);
    ARD = abs(err./he_exp); 
    ind_keep = find(~isinf(ARD));
    ARD= ARD(ind_keep);
    he_pred = he_pred(ind_keep);
    AARD_sys = sum(ARD)/size(he_pred,1);
    AARD_overall(count) = sum(sum(ARD))/prod(size(he_pred));
    filenamesave = strcat('UNIFACerrors',num2str(T),'.mat'); 
    save(filenamesave)
end 
AARD = sum(AARD_overall.*nopreds)/sum(nopreds);
smse = sum(smse_overall.*nopreds)/sum(nopreds);

