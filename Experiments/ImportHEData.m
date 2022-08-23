%% Importing excess enthalpy data into a tensor for use

%Francesca Middleton, 2022-03-02

%% Import the data of composition, component, temperature, real enthalpy and excess enthalpy
clc
clear
data = readtable('HEData9June.xlsx','ReadVariableNames',true); % change sheet to include certain functional groups as the main site of data collection 

comp = table2array(data(:,7));
temp = table2array(data(:,6));
HE  = table2array(data(:,9));
% Unique temperatures 
[B,it,ib]=unique(temp);
count_t = zeros(length(B),2);
B = sort(B, 'descend');
count_t(:,1)=B;
for i = 1:length(B)
    count_t(i,2) = sum(ib(:)==i); 
end 

%% Interpolate data 
% at one temp (ignore pressure) for each mixture to find values at certain concentrations 
% Find all data points with the same temp, FuncGroup1, FuncGroup2,
% ChainLength1 and ChainLength2 and populate comp_temp with compositions
% and HE_temp with HE data 

interp_index = zeros(length(data.FunctionalGroup1),1);%variable to save the indexes that have been interpolated

% Specify the mixtures wanted in the matrix. The algorithm will find all
% combinations of functional group 1 and 2.  
func_groups.one = {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone2','Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc'};
func_groups.two = {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone2','Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc'};
max_chain_length = 12; 
T = 253.15; % temperature in kelvin used for the matrix
% a moderate allowance for different experimental values 
Tupper = T+1;
Tlower = T-1;
P = 15000; % pressure in kPa 
% pressure is ignored due to very small variation with pressure of HE and
% non-critical behaviour 

conc_interval = 0.1:0.1:0.9;
%First restrict data to relevant temperatures. Can add a for loop here
%later 
P_loc = find(data.Pressure_kPa_<P);
temp_data=data(:,:); % all data, do not consider pressure
T_loc = find(temp_data.Temperature<Tupper & temp_data.Temperature>Tlower);
temp_data = temp_data(T_loc, :);

% create matrices to populate with each mixture's data 
HE_data = nan(length(conc_interval), (100*100));
dim1= 216; %assumed max 40 data points per set 
conc_original = nan(dim1, 100*100);
HE_original = nan(dim1, 100*100);
mixture = nan(4, 100*100);
orderPolyfit = nan( 100*100,1);
f1=0;
ind=0;
%run through functional group 1 possibilities 
for func1= func_groups.one
    f1=f1+1;
    f2=0;
    % restrict to mixtures containing the wanted functional group in
    % position 1
    func1_loc= find(strcmp(temp_data.FunctionalGroup1, func1));
    temp1 = temp_data(func1_loc,:);
    % functional group 2 possibilities 
    for func2= func_groups.two
        f2=f2+1;
        % restrict to mixtures containing the wanted functional group in
        % position 2
        func2_loc = find(strcmp(temp1.FunctionalGroup2, func2));
        temp2 = temp1(func2_loc,:);
        % now loop through chain lengths
        for i =0:max_chain_length
            %only this chain length 
            chain1_loc = find(temp2.Chainlength1==i);
            temp3 = temp2(chain1_loc,:);
            for j =0:max_chain_length
                % 2nd chain length considered
                if i~=j || f1~=f2
                    %excess enthalpy of a mixture of the two different pure
                    %components
                    chain2_loc = find(temp3.Chainlength2==j);
                    temp4 = temp3(chain2_loc,:); % this is the data we interpolate 
                    if size(temp4,1)>3
                        ind =ind+1;% this index is correct now 
                        mixture(1,ind)= f1;
                        mixture(2,ind)= i;
                        mixture(3,ind)= f2;
                        mixture(4,ind)= j;
                        % interpolate data and populate matrix of all data
                        [HE_data(:,ind), uncertainty(:,ind), orderPolyfit(ind)]=interp_data(temp4, conc_interval);
                        % save original data in the same order
                        conc_original(1:size(temp4,1),ind) = temp4.Compositioncomponent1;
                        HE_original(1:size(temp4,1),ind) = temp4.Excessenthalpy;
                        % find if mixture is unique
                    end 
                end 
            end  
        end 
    end 
end
 
%remove nan columns or rows 
mixture = mixture(:, all(~isnan(mixture)));
HE_data = HE_data(:, all(~isnan(HE_data)));
conc_original = conc_original(:, all(~isnan(conc_original)));
orderPolyfit = orderPolyfit(all(~isnan(orderPolyfit)), :);
% number of unique components and their indices
[B1, ia1, ib1] = unique(mixture([1,2],:)', 'rows');
[B2, ia2, ib2] = unique(mixture([3,4],:)', 'rows');
% mixture =A: B = A(ia,:) and A = B(ib,:)

%% Populate matrices of interpolated data 
% populate the matrix correctly 
dim1 = size(B1,1); %number of component ones (rows)
dim2 = size(B2,1);% number of component twos (columns)
dim3 = length(conc_interval); % number of compositions

%B1 and B2 already ordered - they become the indices on the sides of the
%matrix 
%order the data in HE_data 
HE_data_sparse = nan(dim1,dim2,dim3);
for i = 1:dim1
    for j=1:dim2
        [lia,index]=ismember([B1(i,:) B2(j,:)], mixture', 'rows');
        if lia ==1
            HE_data_sparse(i,j,:)= reshape(HE_data(:,index),[1,dim3]);
        end 
    end 
end 

%find missing rows and columns 
missing.ind = find(isnan(HE_data_sparse));
% check for nan rows and column 
[missing.i, missing.j] = find(isnan(HE_data_sparse));
% remove all data that is nan
%% Export 3-way array to excel spreadsheet 
filename = strcat('HEMatrix13June',num2str(T),'.xlsx');
%create table with all this information
for i = 1:length(conc_interval)
    Table = array2table(HE_data_sparse(:,:,i));
    writetable(Table,filename,'Sheet',num2str(conc_interval(i)))
end 
Table2= array2table(mixture);
writetable(Table2,filename,'Sheet','mixtures')
TableB1 = array2table(B1);
TableB2 = array2table(B2);
writetable(TableB1,filename, 'Sheet', 'B1')
writetable(TableB2,filename, 'Sheet', 'B2')

%% Check the interpolation worked nicely for select mixtures
clf
mix =58;
plot(conc_original(:,mix), HE_original(:,mix), 'bo')
hold on 
plot(conc_interval, HE_data(:,mix), 'LineWidth', 1) 

legend('Experimental', 'Interpolated', 'Location', 'northwest')

xlabel('Composition component 1 (mol/mol)')
ylabel('Excess enthalpy (kJ/mol)')
title(['Order of the fit ', num2str(orderPolyfit(mix))])
%% Plt the missing data structure 
% these must have the same sizes as x
v=HE_data_sparse;

xslice = 1:1:44;    % location of y-z planes
yslice = 1:2:80;     % location of x-z plane
zslice = 1:3:9;         % location of x-y planes
clf
slice(v,xslice,yslice,zslice)
xlabel('Component 1')
ylabel('Component 2')
zlabel('Composition')
% hm = HeatMap(HE_matrix);
% addXLabel(hm,'Component 1','FontSize',12);
% addYLabel(hm,'Component 2','FontSize',12);
% addZLabel(hm, 'Composition', 'FontSize',12);
% view(hm)
% histogram(X)
% xlabel('Value in the matrix X')
% ylabel('Frequency')


%%
function [HE, uncertainty, orderPolyfit]=interp_data(data, conc_interval)
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
    
    HE_new= HE_original;
    comp_new = comp;
    %interpolate the data 
    maxOrder = 7;
    error = zeros(maxOrder-1,1);
    ind = 1;
    for i =2:maxOrder
        [p,S,mu] = polyfit(comp_new, HE_new,i);
        [~,uncertainty] = polyval(p,conc_interval, S,mu);
        error(ind) = sum(uncertainty)/length(uncertainty);
        ind = ind+1;
    end 
    %populate the data 
    orderPolyfit = find(error == min(error))+1;
    [p,S] = polyfit(comp_new, HE_new,orderPolyfit);
    [HE,uncertainty] = polyval(p,conc_interval, S);

end 
