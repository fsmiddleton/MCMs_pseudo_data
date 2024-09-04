%% Import all collected data 
% Import the data of composition, component, temperature, and excess enthalpy
clc
clear
data = readtable('HEData2023.xlsx','Sheet', 'All','ReadVariableNames',true); % change sheet to include certain functional groups as the main site of data collection 

comp = table2array(data(:,7));
temp = table2array(data(:,6));
HE  = table2array(data(:,9));
rHE = table2array(data(:,11)); %reduced excess enthalpy rHE = HE/(x1*(1-x1))
% Unique temperatures 
[B,it,ib]=unique(temp);
count_t = zeros(length(B),2);
B = sort(B, 'descend');
count_t(:,1)=B;
for i = 1:length(B)
    count_t(i,2) = sum(ib(:)==i); 
end 


f1 = table2cell(data(:,1));
f2 = table2cell(data(:,2));
ch1 = table2array(data(:,3));
ch2 = table2array(data(:,4));

func_groups.one = {'Alkane', 'Primaryalcohol'}%, 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
func_groups.two = {'Alkane', 'Primaryalcohol'}%, 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
max_chain_length = 12; 
f1_num= zeros(length(HE),1);
f2_num= zeros(length(HE),1);
f = 0;
for func = func_groups.one
    f=f+1;
    indices1 = find(strcmp(func,f1));
    indices2 = find(strcmp(func,f2));
    f1_num(indices1) = f;
    f2_num(indices2) = f;
end 
%specification of mixtures 
mixture = zeros(4,length(HE));
mixture(1,:)=f1_num;
mixture(3,:)=f2_num;
mixture(2,:)=ch1;
mixture(4,:)=ch2;
mixtureT = mixture';
%the unique components in the array
[allcomps1,~,~]=unique(mixtureT(:,[1,2]), 'rows');
[allcomps2,~,~]=unique(mixtureT(:,[3,4]), 'rows');
[l,Locb] = ismember(allcomps2,allcomps1,'rows');
include2 = find(Locb==0);
% all possible components in this matrix 
allcomps = [allcomps1; allcomps2(include2,:)];
%% Interpolate data 
% at one temp (ignore pressure) for each mixture to find values at certain concentrations 
% Find all data points with the same temp, FuncGroup1, FuncGroup2,
% ChainLength1 and ChainLength2 and populate comp_temp with compositions
% and HE_temp with HE data 
poly = 1; %1 for interpolation using polynomials, 2 for the RK equation, done with reduced HE data 
interp_index = zeros(length(data.FunctionalGroup1),1);%variable to save the indexes that have been interpolated


max_chain_length = 12; 
P = 15000; % pressure in kPa 
% pressure is ignored due to very small variation with pressure of HE and
% non-critical behaviour 

conc_interval = 0.05:0.05:0.95;
Temps = [288.15; 298.15; 313.15; 318.15;323.15]; %[288.15; 298.15; 303.15; 307.5; 309.5;313.15;318.15; 323.15]; % ;243.15; 253.15; 263.15; 273.15; 283.15; 288.15; 290.15; 293.15; 296.15; 298.15; 303.15; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15; 363.15];
Tind = 0;


comp = table2array(data(:,7));
temp = table2array(data(:,6));
HE  = table2array(data(:,9));
rHE = table2array(data(:,11));



%First restrict data to relevant temperatures. Can add a for loop here
%later 
P_loc = find(data.Pressure_kPa_<P);

dimarray = 1000;
dim1= 300; %assumed 
dim2 = dim1;
dim3 = length(conc_interval);
dim4 = length(Temps);
% create matrices to populate with each mixture's data 
HE_data = nan(length(conc_interval), dimarray, length(Temps));
HEpred = nan(dim1,  dimarray,length(Temps));
errorpred = nan(dim1,  dimarray,length(Temps));
uncertainty = nan(length(conc_interval), dimarray,length(Temps));

conc_original = nan(dim1,  dimarray,length(Temps));
HE_original = nan(dim1,  dimarray,length(Temps));
mixture = nan(4,  dimarray);
mxture2 = nan(4,dimarray);
TempsMix = zeros(length(Temps), dimarray);
orderPolyfit = nan(  dimarray,1);
f1=0;
ind=0;
disp('Regressing')
%4-way array



% 3-way array 
%run through functional group 1 possibilities 
for func1= func_groups.one
    f1=f1+1;
    f2=0;
    % restrict to mixtures containing the wanted functional group in
    % position 1
    func1_loc= find(strcmp(data.FunctionalGroup1, func1));
    temp1 = data(func1_loc,:);
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

                    if size(temp4,1)>2 
                        %restricted to mixture wanted 
                        ind =ind+1;% this index is correct now 

                        mixture(1,ind)= f1;
                        mixture(2,ind)= i;
                        mixture(3,ind)= f2;
                        mixture(4,ind)= j;
                        % interpolate data and populate matrix of all data
                        Tind = 0;
                        for T = Temps'
                            Tind = Tind+1;
                            % a moderate allowance for different experimental values 
                            Tupper = T+1;
                            Tlower = T-1;
                            T_loc = find(temp4.Temperature<Tupper & temp4.Temperature>Tlower);
                            temp5 = temp4(T_loc, :);
                            if size(temp5,1)>2
                                TempsMix(Tind,ind) = 1; % a 1 is shown in the array for every existing temperature of a mixture 
                                if poly == 1
                                    [HE_data(:,ind, Tind), uncertainty(:,ind,Tind), orderPolyfit(ind,Tind),HEpred(:,ind,Tind), errorpred(:,ind,Tind)]=interp_data(temp5, conc_interval);
                                    HE_original(1:size(temp5,1),ind,Tind) = temp5.Excessenthalpy;
                                else
                                    %used for reduced HE data 
                                    [HE_data(:,ind, Tind), uncertainty(:,ind,Tind), orderPolyfit(ind,Tind),HEpred(:,ind,Tind), errorpred(:,ind,Tind)]=RK_interp_data(temp5, conc_interval);
                                    HE_original(1:size(temp5,1),ind, Tind) = temp5.ReducedHE;
                                end 
                                % save original data in the same order
                                conc_original(1:size(temp5,1),ind, Tind) = temp5.Compositioncomponent1;
                            end 
                        end 
                    end 
                end  
            end 
        end 
    end 
end  


% Create 4-way array 
    disp('Exporting')
    prefixfilename = 'HE4wayArrayPolySmall5';
    %remove nan columns or rows 
    mixture = mixture(:, 1:ind);
    HE_data = HE_data(:, 1:ind,:);
    conc_original = conc_original(:, 1:ind,:);
    HE_original = HE_original(:, 1:ind,:);
    uncertainty = uncertainty(:, 1:ind,:);
    %orderPolyfit = orderPolyfit(1:ind, :)';
    RAD = abs(errorpred./HEpred);
    % number of unique components and their indices
    mixtureT = mixture';
    [comps1,~,~]=unique(mixtureT(:,[1,2]), 'rows');
    [comps2,~,~]=unique(mixtureT(:,[3,4]), 'rows');
    [l,Locb] = ismember(comps2,comps1,'rows');
    include2 = find(Locb==0);
    % all possible components in this matrix 
    comps = [comps1; comps2(include2,:)];

    % Populate matrices of interpolated data 
    % populate the square array 
    dim1 = size(comps,1); %number of component ones (rows)
    dim2 = size(comps,1);% number of component twos (columns)
    dim3 = length(conc_interval); % number of compositions
    dim4 = length(Temps);

    %B1 and B2 already ordered - they become the indices on the sides of the
    %matrix 
    %order the data in HE_data 
    HE_data_sparse = nan(dim1,dim2,dim3,dim4);
    error = nan(dim1,dim2,dim3,dim4);
    for i = 1:dim1
        for j=1:dim2
            for k =1:dim3
                %Find mixture
                [lia,index]=ismember([comps(i,:) comps(j,:)], mixture', 'rows');
                if lia ==1
                    %Find temps 
                    Temperaturestemp = TempsMix(:,index);
                    HE_data_sparse(i,j,:,:)= reshape(HE_data(:,index,:),[1,dim3,dim4]);
                    HE_data_sparse(j,i,:,:)=reshape(flip(HE_data(:,index,:)),[1,dim3,dim4]);
                    error(i,j,:,:) = reshape(uncertainty(:,index,:),[1,dim3,dim4]);
                    error(j,i,:,:) = reshape(uncertainty(:,index,:),[1,dim3,dim4]);
                    
                end 
            end 
        end 
    end 
    %place zeros on the diagonal 
    for i =1:dim1
        HE_data_sparse(i,i,:,:)=zeros(length(conc_interval),length(Temps));
    end 
    %find missing rows and columns 
    missing.ind = find(isnan(HE_data_sparse));
    % check for nan rows and column 
    [missing.i, missing.j] = find(isnan(HE_data_sparse));
    observed = find(~isnan(HE_data_sparse));
    total = length(missing.ind)+length(observed);
    totalpercobserved = length(observed)/total;
    observed3 = zeros(length(Temps),1);
    total3 = dim1*dim2*dim3;
    for indexobs = 1:length(Temps)
        observed3(indexobs) = length(find(~isnan(HE_data_sparse(:,:,:,indexobs))));
        percobs3(indexobs) = observed3(indexobs)/total3*100;
    end 
    filename = strcat(prefixfilename,'.mat');
    save(filename)
    disp('Exported')
%% Save array 
save HEData4waySmallPoly 

%% Export 4-way array 
    % Export 4-way array to excel spreadsheet 
    filename = strcat(prefixfilename,'.xlsx');
    %create table with all this information
    for T = 1:length(Temps)
        for i = 1:length(conc_interval)
            Table = array2table(HE_data_sparse(:,:,i));
            writetable(Table,filename,'Sheet',strcat(num2str(Temps(T)),'-',num2str(conc_interval(i))))
        end 
    end 
    TableHE = array2table(HE_data);
    Table2= array2table(mixture);
    Table4 = array2table(uncertainty);
    Table5 = array2table(orderPolyfit);
    TableRAD = array2table(RAD);
    TableConcO = array2table(conc_original);
    TableHEO = array2table(HE_original);
    writetable(Table2,filename,'Sheet','mixtures1')
    writetable(Table4, filename, 'Sheet', 'error')
    writetable(Table5, filename, 'Sheet', 'orderPolynomial')
    TableB1 = array2table(comps1); %unique component 1
    TableB2 = array2table(comps2);%unique component 2
    writetable(TableB1,filename, 'Sheet', 'B1')
    writetable(TableB2,filename, 'Sheet', 'B2')
    writetable(TableRAD, filename, 'Sheet', 'RAD')
    writetable(TableConcO, filename, 'Sheet', 'ConcOriginal')
    writetable(TableHEO, filename, 'Sheet', 'HEOriginal')
    writetable(TableHE, filename, 'Sheet', 'HEInterpolated')
    disp(filename)
    disp('Exported')

%%
function [HE, uncertainty, orderPolyfit, HEpredout, errorpredout]=interp_data(data, conc_interval)
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
    % check for data that does not have many points outside a certain interval 
    maxdiff = max(comp)-min(comp);
    numberofpoints = length(comp_new);
    if ind_keep
        comp_new(numberofpoints+1)=0;
        comp_new(numberofpoints+2)=1;
        HE_new(numberofpoints+1)=0;
        HE_new(numberofpoints+2)=0;
    end 
    
    %interpolate the data 
    if maxdiff <0.3 || numberofpoints <2
        %just fit a parabola to the small amount of data 
        maxOrder = 2;
    else 
        maxOrder = 7;
    end 
    
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
    [HEpred, errorpred] = polyval(p, comp_new,S);
    HEpredout = zeros(300,1);
    HEpredout(1:length(HEpred),1)=HEpred;
    errorpredout = zeros(300,1);
    errorpredout(1:length(errorpred),1)=errorpred;
end 

function [HEr, uncertainty, orderPolyfit, HEpredout, errorpredout]=RK_interp_data(data, conc_interval)
    % Interpolate data for use in the matrices using the Redlich Kister
    % equation:
    % HE /(x1(1-x))=( A+B(2x-1)+C(2x-1)^2+D(2x-1)^3+E(2x-1)^4+F(2x-1)^5)
    % Inputs specify the mixture to be captured and the concentration
    % interval for which to create interpolate points
    
    % data = the whole table of data to be interpolated
    % conc_interval = list of concentrations for interpolation eg conc_interval = 0.1:0.1:0.9;
    
    % Ouputs is the data produced by the interpolation, ready for
    % matrix populating 
    % HE = reduced excess enthalpy data at each concentration value for the mixture
    % uncertainty = uncertainty associated with each data point, to be
    % added
    % order = the number of constants used in the RK equation
    HEr_original = data.ReducedHE;
    comp = data.Compositioncomponent1;
    % remove the 0 and 1 data points 
    ind_keep = find(HEr_original ~=0);
    HEr_new = HEr_original(ind_keep);
    comp_new = comp(ind_keep);
    maxdiff = max(comp)-min(comp);
    numberofpoints = length(comp_new);
    if ind_keep
        comp_new(numberofpoints+1)=0;
        comp_new(numberofpoints+2)=1;
        HEr_new(numberofpoints+1)=0;
        HEr_new(numberofpoints+2)=0;
    end 
    
     %interpolate the data 
    if maxdiff <0.4 || numberofpoints<4
        %just fit a parabola to the small amount of data 
        maxOrder = 3;
    else 
        maxOrder = 7;
    end 
    %interpolate the data 
    
    error = zeros(maxOrder-1,1);
    ind = 1;
    for i =2:maxOrder
        %Y=XB, B=X\Y, Y=HEr, X= (2x-1)^i, B=constants
        Y=HEr_new; % column n-by-1
        X = ones( length(comp_new),i); %
        for j = 1:i
            X(:,j)=(2*comp_new-1).^(j-1);
        end 
        
        mdl = fitlm(X,Y, 'Intercept',false);
        coeff = mdl.Coefficients.Estimate;
       
        error(ind) = mdl.MSE;
        ind = ind+1;
    end 
    %Find the best fit 
    orderPolyfit = find(error == min(error))+1;
    if length(orderPolyfit)>1
        orderPolyfit = orderPolyfit(1,1);
    end
    X = ones(length(comp_new), orderPolyfit); %
    for j = 1:i
        X(:,j)=(2*comp_new-1).^(j-1);
    end 
    mdl = fitlm(X,Y,'Intercept',false);
    coeff = mdl.Coefficients.Estimate;
    %predict the data 
    Xpred = ones(length(conc_interval), orderPolyfit); %
    for j = 1:i
        Xpred(:,j)=(2*conc_interval-1).^(j-1);
    end 
    HEpred = ones(length(conc_interval), orderPolyfit); %
    for j = 1:i
        HEpred(:,j)=(2*conc_interval-1).^(j-1);
    end  
    HEr = Xpred*coeff;
    HEpred = HEpred*coeff;
    errorpred = HEpred - HEr_new;
    uncertainty = mdl.MSE;
    HEpredout = zeros(300,1);
    HEpredout(1:length(HEpred),1)=HEpred;
    errorpredout = zeros(300,1);
    errorpredout(1:length(errorpred),1)=errorpred;
end 
