% Create 4-way arrays 
% FS Middleton 27 June 2022

%% Import data 
clc
clear
data = readtable('HEData2023.xlsx','Sheet', 'All','ReadVariableNames',true); % change sheet to include certain functional groups as the main site of data collection 

comp = table2array(data(:,7));
temp = table2array(data(:,6));
HE  = table2array(data(:,9));
rHE = table2array(data(:,11)); %reduced excess enthalpy rHE = HE/(x1*(1-x1))...
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

func_groups.one = {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
func_groups.two = {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
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
allcomps = [allcomps1; allcomps2(include2,:)];
%% Interpolate data 

poly = 1; %1 for interpolation using polynomials, 2 for the RK equation, done with reduced HE data 
interp_index = zeros(length(data.FunctionalGroup1),1);%variable to save the indexes that have been interpolated
max_chain_length = 12; 
P = 15000; % pressure in kPa 
conc_interval = 0.05:0.05:0.95;
Temps = [293.15;298.15; 303.15;308.15;313.15;318.15;323.15];
Tind = 0;
comp = table2array(data(:,7));
temp = table2array(data(:,6));
HE  = table2array(data(:,9));
rHE = table2array(data(:,11));
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
%4-way array
for func1= func_groups.one
    f1=f1+1;
    f2=0;
    % restrict to mixtures containing the desired functional group in
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
                        Tind = 0;
                        for T = Temps'
                            Tind = Tind+1;
                            % a moderate allowance for different experimental values 
                            Tupper = T+1;
                            Tlower = T-1;
                            T_loc = find(temp4.Temperature<Tupper & temp4.Temperature>Tlower);
                            temp5 = temp4(T_loc, :);
                            if size(temp5,1)>2
                                ind =ind+1;
                                mixture(1,ind)= f1;
                                mixture(2,ind)= i;
                                mixture(3,ind)= f2;
                                mixture(4,ind)= j;
                                TempsMix(Tind,ind) = 1;
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
prefixfilename = strcat('HE4wayArrayPolySmall');%,num2str(length(Temps)));
%remove nan columns or rows 
mixture = mixture(:, 1:ind);
HE_data = HE_data(:, 1:ind,:);
conc_original = conc_original(:, 1:ind,:);
HE_original = HE_original(:, 1:ind,:);
uncertainty = uncertainty(:, 1:ind,:);
RAD = abs(errorpred./HEpred);
comps = allcomps;

% Populate matrices of interpolated data 
% populate the square array 
dim1 = size(comps,1); %number of component ones (rows)
dim2 = size(comps,1);% number of component twos (columns)
dim3 = length(conc_interval); % number of compositions
dim4 = length(Temps);
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
filename = strcat(prefixfilename,'.mat');
save(filename)
