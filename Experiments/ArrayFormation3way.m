% Create 3-way arrays 
% FS Middleton 27 June 2022

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

func_groups.one = {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
func_groups.two = {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
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
poly = 1; %1 for interpolation using polynomials, 2 for the RK equation, done with reduced HE data 
interp_index = zeros(length(data.FunctionalGroup1),1);%variable to save the indexes that have been interpolated
P = 15000; % pressure in kPa 
%% Interpolate data 
for interval = 0.05
    conc_interval = interval:interval:(1-interval);
    Temps = [293.15; 298.15; 303.15; 308.15; 313.15; 318.15; 323.15];
    for T = Temps'
        data = readtable('HEData2023.xlsx','ReadVariableNames',true);
        comp = table2array(data(:,7));
        temp = table2array(data(:,6));
        HE  = table2array(data(:,9));
        rHE = table2array(data(:,11));
        % a moderate allowance for different experimental values 
        Tupper = T+1;
        Tlower = T-1;
        P_loc = find(data.Pressure_kPa_<P);
        temp_data=data(:,:); % all data, do not consider pressure
        T_loc = find(temp_data.Temperature<Tupper & temp_data.Temperature>Tlower);
        temp_data = temp_data(T_loc, :);
        % create matrices to populate with each mixture's data 
        HE_data = nan(length(conc_interval), (100*100));
        HEpred = nan(300, (100*100));
        errorpred = nan(300, (100*100));
        uncertainty = nan(length(conc_interval),100*100);
        dim1= 216; %assumed max data points per set 
        conc_original = nan(dim1, 100*100);
        HE_original = nan(dim1, 100*100);
        mixture = nan(4, 100*100);
        mxture2 = nan(4,100*100);
        orderPolyfit = nan( 100*100,1);
        R2 = nan(100*100);
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
                % restrict to mixtures containing functional group inposition 2
                func2_loc = find(strcmp(temp1.FunctionalGroup2, func2));
                temp2 = temp1(func2_loc,:);
                %loop through chain lengths
                for i =0:max_chain_length
                    chain1_loc = find(temp2.Chainlength1==i);
                    temp3 = temp2(chain1_loc,:);
                    for j =0:max_chain_length
                        if i~=j || f1~=f2 
                            chain2_loc = find(temp3.Chainlength2==j);
                            temp4 = temp3(chain2_loc,:); % this is the data we interpolate 
                            if size(temp4,1)>2
                                ind =ind+1;% this index is correct now
                                mixture(1,ind)= f1;
                                mixture(2,ind)= i;
                                mixture(3,ind)= f2;
                                mixture(4,ind)= j;
                                % interpolate data and populate matrix of all data
                                if poly == 1
                                    [HE_data(:,ind), uncertainty(:,ind), orderPolyfit(ind),HEpred(:,ind), errorpred(:,ind), R2(ind)]=interp_data(temp4, conc_interval);
                                    HE_original(1:size(temp4,1),ind) = temp4.Excessenthalpy;
                                else
                                    %used for reduced HE data 
                                    [HE_data(:,ind), uncertainty(:,ind), orderPolyfit(ind),HEpred(:,ind), errorpred(:,ind)]=RK_interp_data(temp4, conc_interval);
                                    HE_original(1:size(temp4,1),ind) = temp4.ReducedHE;
                                end 
                                % save original data in the same order
                                conc_original(1:size(temp4,1),ind) = temp4.Compositioncomponent1;
                                % find if mixture is unique
                            end 
                        end 
                    end  
                end 
            end 
        end
        % Exporting the data 
        prefixfilename = strcat('HE3wayPoly', num2str(T));
        %remove nan columns or rows 
        mixture = mixture(:, 1:ind);
        HE_data = HE_data(:, 1:ind);
        conc_original = conc_original(:, 1:ind);
        uncertainty = uncertainty(:, 1:ind);
        orderPolyfit = orderPolyfit(1:ind, :)';
        RAD = abs(errorpred./HEpred);
        % number of unique components and their indices
        mixtureT = mixture';
        [comps1,~,~]=unique(mixtureT(:,[1,2]), 'rows');
        [comps2,~,~]=unique(mixtureT(:,[3,4]), 'rows');
        [l,Locb] = ismember(comps2,comps1,'rows');
        include2 = find(Locb==0);
        % all possible components in this matrix 
        comps = union(comps1, comps2, 'rows'); 
        dim1 = size(comps,1); %number of component ones (rows)
        dim2 = size(comps,1);% number of component twos (columns)
        dim3 = length(conc_interval); % number of compositions
        %Create the 3-way array
        HE_data_sparse = nan(dim1,dim2,dim3);
        error = nan(dim1,dim2,dim3);
        count =0;
        for i = 1:dim1
            for j=1:dim2
                %consider all possibilities of mixtures 
                [lia,index]=ismember([comps(i,:) comps(j,:)], mixture', 'rows');
                if lia ==1
                    count = count+1;
                    HE_data_sparse(i,j,:)= reshape(HE_data(:,index),[1,dim3]);
                    HE_data_sparse(j,i,:)=reshape(flip(HE_data(:,index)),[1,dim3]);
                    error(i,j, :) = reshape(uncertainty(:,index),[1,dim3]);
                    error(j,i,:) = reshape(uncertainty(:,index),[1,dim3]);
                    if i==j
                        disp(i)
                    end 
                end 
            end 
        end 
        %place zeros on the diagonal 
        for i =1:dim1
            HE_data_sparse(i,i,:)=zeros(length(conc_interval),1);
        end 
        save(strcat(prefixfilename,'.mat'))
    end 
end 