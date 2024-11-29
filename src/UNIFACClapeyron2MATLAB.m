%% Importing UNIFAC predictions 
% Imports UNIFAC predictions made in Clapeyron to MATLAB and writes into a
% MATLAB file 

%% Import data 
clc
clear
%%
% choose which temperature you wish to transfer to matlab
Temps =[293.15; 298.15; 303.15;308.15; 313.15; 318.15; 323.15];
func_group_no= {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
%import compound data 
compounds = readtable("data/UNIFACParams.xlsx", 'Sheet', "compounds");
functionalgroup = table2cell(compounds(:,1));
chainlength = table2array(compounds(:,2));
compounds = table2array(compounds(:,3));
%compounds{84} = 'tertbutylalcohol';
conc_interval = 0:0.01:1;

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
        mixture(1,i) = find(strcmp(functionalgroup(ind1), func_group_no));
        mixture(3,i) = find(strcmp(functionalgroup(ind2), func_group_no));
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

%initiate necessary vars 
Temps =[313.15]%; 298.15; 303.15;308.15; 313.15; 318.15; 323.15];
smse_overall = zeros(length(Temps),1);
nopreds = zeros(length(Temps),1);
AARD_overall = zeros(length(Temps),1);
counter=0;% temperature counter 
for T = Temps'% for each excel file
    counter=counter+1;
    %import experimental data 
    filename = strcat('HE3wayPolyAll',num2str(T),'.mat');
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
    conc_indices = 6:5:96; %known from previous code - indices of predictions 
    ind_keep = 1:size(he_exp,2);
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
    mixture_exp = mixture_exp(:,find(ind_keep));
    %put the mixtures in the correct order 
    for mix = 1:size(mixture_exp,2)
        if mixture_exp(1,mix)<mixture_exp(3,mix)
            mixture_exp(:,mix) = [mixture_exp(3:4,mix); mixture_exp(1:2,mix)];
        end 
    end 
    err = he_exp-he_pred;
    
    %remove zeros 
    smse_sys = sqrt(sum(err.^2)/size(he_pred,1)); % per system
    smse_overall(counter) = sqrt(sum(sum(err.^2))/prod(size(he_pred)));
    nopreds(counter) = size(he_pred,1);
    ARD = abs(err./he_exp); 
    ind_keep = find(~isinf(ARD));
    ARD = reshape(ARD(ind_keep), [19], []);
    err = reshape(err(ind_keep), [19], []);
    he_expard = reshape(he_exp(ind_keep), [19], []);
    he_predard = reshape(he_pred(ind_keep), [19],[] );
    AARD_sys = sum(ARD)/size(he_predard,1);
    AARD_overall(counter) = sum(sum(ARD))/prod(size(he_predard));
    filenamesave = strcat('UNIFACerrors',num2str(T),'.mat'); 
    save(filenamesave)
end 
AARD = sum(AARD_overall.*nopreds)/sum(nopreds);
smse = sum(smse_overall.*nopreds)/sum(nopreds);

%% Parity plot and error distribution plots 
%Import all the UNIFAC predictions 
Temperature = 298.15;%[288.15; 298.15; 303.15;307.5;309.5; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15];
load('colorblind_colormap.mat')

load(strcat('UNIFACerrors',num2str(Temperature),'.mat')) %change the temperature here
 
%parity plot 
f4 = figure(4);
f4.Position=[10 10 350 250];
clf
%choose between plot and loglog 
hold on
plot(he_exp(:), he_pred(:),'.','Color', 'k', 'MarkerSize',6)
plotx = he_exp;
plotx(isnan(plotx)) =[];
%prediction interval
PIx = min(plotx(:)):20:max(plotx(:));
syx = std(plotx(:));
hi = 1/numel(plotx(:)) + (PIx-mean(PIx)).^2/sum((PIx-mean(PIx)).^2);
PIy1 = PIx - 1.96.*syx.*sqrt(1+hi);
PIy2 = PIx + 1.96.*syx.*sqrt(1+hi);
plot(PIx,PIy1, '--', 'Color', colorblind(7,:), 'LineWidth', 1) %make this pink 
plot(PIx,PIy2, '--', 'Color', colorblind(7,:), 'LineWidth', 1) %make this pink 
%y=x line 
plot([min(plotx(:)),max(plotx(:))], [min(plotx(:)),max(plotx(:))],'-', 'Color',colorblind(4,:), 'LineWidth', 1.2)
legend('UNIFAC','PI (95%)', 'FontSize',8,'Location', 'northwest');
hold off
xlabel('Experimental (J/mol)','FontSize',11, 'Color', 'k')
ylabel('UNIFAC prediction (J/mol)','FontSize',11, 'Color', 'k')

%% error distribution 
f5 = figure(5);
f5.Position = [10 20 350 250];
clf
%choose the number of bins 
%ploterrhist(errorHE(:), 'bins', 20)
histogram(err(:), 20, 'FaceColor', colorblind(6,:), 'BinLimits', [-1000 1000])
ylabel('Instances','FontSize',11, 'Color', 'k')
xlabel('UNIFAC error (J/mol)', 'FontSize',11, 'Color', 'k')

%% 3way error distribution 
f6 = figure(6);
f6.Position = [10 20 350 250];
clf
edges = {-1000:100:1000; 0:5:100}; % bin edges
data = [err(:) abs(ARD(:)).*100]; 
hist3(data, 'Edges', edges,'CdataMode','auto','FaceColor','interp')
zlabel('Instances','FontSize',11, 'Color', 'k')
xlabel('Error (J/mol)', 'FontSize',11, 'Color', 'k')
ylabel('ARD (%)' ,'FontSize', 11, 'Color', 'k')
%colormap color
colormap(colorblind(flip([9,1,8,6,3]),:))
%colormap(flip(parula))
% Exporting figure 
%saveas(f6, 'test.png');
set(f6, 'Color', 'w');
export_fig UNIFACHist3-318.15.jpg -nocrop -native
% ref: Yair Altman (2023). export_fig (https://github.com/altmany/export_fig/releases/tag/v3.29), GitHub. Retrieved January 12, 2023.
%% Errors by type of mixture 
clc
clear

he_predtype = cell(14,5);
he_exptype = cell(14,5);
count = zeros(14,5);
counter= 0; % temperature counter;%
for T =[293.15; 298.15; 303.15; 308.15; 313.15; 318.15; 323.15]'
    counter = counter+1;
    filename = strcat('UNIFACerrors',num2str(T),'.mat');
    load(filename)
    % define types of mixtures 
    dim = size(he_exp);
    type = zeros(1,dim(2));
    for i = 1:dim(2)
        %extract data set 

        temp = he_exp(:,i);
        maxVal = max(temp);
        minVal = min(temp);
        if maxVal<0 || minVal<0 
            % passes through the origin or is negative 
            if sign(maxVal)~=sign(minVal) %passes through origin, only one is negative 
                type(i) = 1;
            else 
                type(i) = 2;
            end 
        else 
            %positive curve 
            if maxVal > 1000 % adjustable ceilings for this classification
                type(i) = 3; % very positive 
            elseif maxVal >200
                type(i) = 4; % moderately positive 
            else 
                type(i) = 5; % less than 200 at the maximum 
            end 
        end 
    end
    for i = 1:5
        he_predtype{counter,i} = he_pred(:,find(type==i));
        he_exptype{counter,i} = he_exp(:,find(type==i));
        countofType(counter,i) = length(find(type==i));
    end 
    save(filename)
end 
save('UNIFACPredsPerType.mat')

%%
load('UNIFACPredsPerType.mat')
he_exppertype = cell(5,1);
he_predpertype = cell(5,1);
for i =1:5
    he_exps = [];
    he_preds = [];
    for j = 1:7 %number of temperatures 
        he_exps = [he_exps he_exptype{j,i}];
        he_preds = [he_preds he_predtype{j,i}];
    end 
    he_exppertype{i,1} = he_exps;
    he_predpertype{i,1} = he_preds;
end
save('UNIFACPredsPerType.mat')
%% Tables per type of mixture 
load('UNIFACPredsPerType.mat')
Temps = [293.15; 298.15; 303.15; 308.15; 313.15; 318.15; 323.15]';
%choose the temperature 
T = 293.15; 
indexT = find(Temps==T);
for i =1:5 
    numberofMix(i) = countofType(indexT,i);
    heexp = he_exptype{indexT,i};
    hepred = he_predtype{indexT,i};
    hepred(heexp==0)=[];
    heexp(heexp==0) = [];
    
    err = heexp(:)-hepred(:);
    err(isnan(err))=[];
    heexp = heexp(:);
    heexp(isnan(heexp))=[];
    
    mse(i) = sqrt(sum(err(:).^2)/length(err(:)));
    wmse(i) = sqrt(find_wmse_error(err(:), length(err(:))));
    ard(i) = mean(abs(err(:)./heexp(:)))*100;
end 

tbl = table(["1","2","3","4","5"]',numberofMix',ard',mse',wmse', 'VariableNames', ["Type of mixture", "Number of mixtures","AARD (%)","SMSE (J/mol)","wSMSE (J/mol)" ]);
summary(tbl)

%% Errors per type
load('UNIFACPredsPerType.mat')
error_type = cell(5,1);
ARD_type = cell(5,1);
ARD_typeoverall = zeros(5,1);
smse_typeoverall = zeros(5,1);
for i = 1:5
    preds = he_predpertype{i,1};
    exps = he_exppertype{i,1};
    err = preds-exps;
    ARD = abs(err./exps);
    ARD = ARD(find(~isinf(ARD)));
    error_type{i} = err;
    ARD_type{i} = ARD;
    ARD_typeoverall(i) = sum(ARD(:))/length(ARD(:));
    smse_typeoverall(i) = sqrt(sum(err(:).^2)/length(err(:)));
end 
save('UNIFACPredsPerType.mat')


%% Errors by functional groups 
clc
clear
Temps =298.15%[288.15; 298.15; 303.15;307.5;309.5; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15];
%Temps = 298.15;%[288.15; 298.15; 303.15;307.5;309.5; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15];
%func_group= {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
func_group_no = {1,2,3:5, [6,7,11], 8:9, 10, 12:17,18:20, 21, 22};
label_func = {'Alkane', 'Primary alcohol', 'Other alcohol','Cycloalkanes', 'Ketone', 'Alkene', 'Ester', 'Amine','Acid','Aldehyde'};
%import preds 

err_funcgroup = cell(length(label_func ),length(label_func ),length(Temps));
ARD_funcgroup = cell(length(label_func ),length(label_func ),length(Temps));


for T = Temps'
    filename = strcat('UNIFACerrors',num2str(T),'.mat');
    load(filename)
    for funcind1 = 1:length(label_func)
        for funcind2 = 1:length(label_func)
            %extract arrays of groups 
            funcgroup1 = func_group_no{funcind1};
            funcgroup2 = func_group_no{funcind1};
            %find mixtures with this combination
            index = find(mixture_exp(1,:)==funcind1);
            index2 = find(mixture_exp(3,:)==funcind2);
            indices = intersect(index,index2);
            %export errors and ARDs
            errtemp = err(:,indices);
            err_funcgroup{funcind1,funcind2} = errtemp; 
            smse_funcgroup(funcind1,funcind2) = sqrt(sum(errtemp(:).^2)/length(errtemp(:)));
            wsmse_funcgroup(funcind1,funcind2) = find_wmse_error(errtemp(:),length(errtemp(:)));
            ardtemp = ARD(:,indices);
            ARD_funcgroup{funcind1,funcind2} = ardtemp;
            AARD_funcgroup(funcind1,funcind2) = sum(ardtemp(:))/length(ardtemp(:)); 
            nomix(funcind1,funcind2) = size(errtemp,2);
        end 
    end 
    save(filename)
end 

%% heatmap of errors by functional group 
clc 
clear
T=298.15;
filename = strcat('heUNIFACforT=',num2str(T),'.mat');
load(filename)
[X,Y] = meshgrid(1:size(AARD_funcgroup,1), 1:size(AARD_funcgroup,1));
SMSEplot = smse_funcgroup;
AARDplot=AARD_funcgroup.*100; 
%create a table 
[row,col] = find(~isnan(AARDplot));
for i = 1:length(row)
    funcgroup1tbl{i} = label_func{row(i)};
    funcgroup2tbl{i} = label_func{col(i)};
    smse_tbl(i) = SMSEplot(row(i), col(i));
    aard_tbl(i) = AARDplot(row(i),col(i));
    wsmse_tbl(i) = wsmse_funcgroup(row(i),col(i));
end 
tbl = table(funcgroup1tbl', funcgroup2tbl',wsmse_tbl',smse_tbl',aard_tbl', 'VariableNames', ["Functional group 1", "Functional group 2", "wSMSE (J/mol)", "SMSE (J/mol)","AARD (%)" ]);
% heatmap 
f1 = figure(1);
f1.Position=[10 10 1000 400];
clf
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', "AARD (%)" , 'ColorMethod', 'none');
h.Colormap=winter;
h.FontSize = 14;

f2 = figure(2);
f2.Position=[10 10 1000 400];
clf
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'SMSE (J/mol)', 'ColorMethod', 'none');
h.Colormap=winter;
h.FontSize = 14;

f3 = figure(3);
f3.Position=[10 10 1000 400];
clf
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'wSMSE (J/mol)', 'ColorMethod', 'none');
h.Colormap=winter;
h.FontSize = 14;
