%% Francesca Middleton 27/09/2022

%discovering the percent correct prediction of the sign predictions 

clc
clear
load('colorblind_colormap.mat')
load('3wayHOSVD-threshold60-fn3=1-Small-concslices-LOOCV-Xnone-maxiter=20000-Temp=318.15-fillmethod=row-21-Dec-2022.mat')
%load('3wayHOSVD-threshold60-centre+scale-Small-concslices-LOOCV-Xnone-maxiter=20000-Temp=318.15-fillmethod=avr-09-Dec-2022.mat')
%% No sign predictions exist
load("2waySVDUnfoldingMode1-21comps-threshold60par-LOOCV-Xscale-T=298.15-fillmethod=reg-21-Nov-2022.mat")
%load('3wayHOSVDThreshold-centre+scale-Fillchanged-A1-Small-concslices-LOOCV-Xscale-maxiter=20000-Temp=298.15-fillmethod=reg-17-Nov-2022.mat')
%load('3wayHOSVDThreshold-Nocentre+scale-Fillchanged-A3-Small-concslices-LOOCV-Xscale-maxiter=20000-Temp=298.15-fillmethod=reg-17-Nov-2022.mat')
%load('3wayHOSVDThresholdsmallercentre+scaleFillchangedA3-Small-concslices-LOOCV-Xscale-maxiter=20000-Temp=298.15-fillmethod=reg-16-Nov-2022.mat')
%load('3wayHOSVDThresholdCentre+ScaleA1-Small-concslices-LOOCV-Xscale-maxiter=20000-Temp=298.15-fillmethod=avr-16-Nov-2022.mat')
%load('3wayHOSVDThresholdCentreA1-Small-concslices-LOOCV-Xscale-maxiter=20000-Temp=298.15-fillmethod=avr-15-Nov-2022.mat');

%Options
Temperature=298.15;
%if only one side of the 2-way array is considered or both (1:2)
whichside =1:2;
postprocess = 0;
Tempind =2;


load(filename)
mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 

%finds the true sign values
X(X==0)=nan;
% for counter =1:size(X,3)
%     Xnew(:,:,counter,Tempind) = remove_nan2(X(:,:,counter,Tempind));
% end 
%X = Xnew;
Xsign = sign(X);
Xscale = Xsign.*log(Xsign.*X);
concentrations = conc_interval;
tempX = tril(X(:,:,1,Tempind),-1)+triu(nan(size(X(:,:,1,Tempind))));
[row,col] = find(~isnan(tempX));
filled_ind =  find(~isnan(tempX));
         
clear Truthsign
for index = 1:length(filled_ind)
    Truthsign(:,index,1) = Xsign(row(index),col(index),:,Tempind);
    Truthsign(:,index,2) = (Xsign(col(index),row(index),:,Tempind)); 
end 



%% Xscale
date ='13-Oct';

scalePrediction = Xm_boot;
scalePrediction(scalePrediction ==0)=nan;
indicesScale = find(~isnan(scalePrediction(1,1,:)));
%indicesInt = intersect(indicesSign,indicesScale); % k that you can use to judge He predictions 
Xtemp = Xs(:,:,1,Tempind);
tempX = tril(Xtemp,-1)+triu(nan(size(Xtemp)));
clear Truthscale

X(X==0)=nan;
% for counter =1:size(X,3)
%     Xnew(:,:,counter) = remove_nan2(X(:,:,counter));
% end 

%X = Xnew;
Xscale = log(sign(X).*(X)).*sign(X);
scalePreds = zeros(length(fns),size(X,3),length(filled_ind),2);
mse = zeros(length(fns),size(X,3));
wmse = zeros(length(fns),size(X,3));
aard = zeros(length(fns),size(X,3));
ARD = zeros(length(fns),size(X,3),length(filled_ind),2);
errors = zeros(length(fns),size(X,3),length(filled_ind),2);

Preds = zeros(length(concentrations),length(filled_ind),length(whichside));
Xm_boot = reshape(Xm_boot,[],length(fns),length(filled_ind));
Xm_boot2 = reshape( Xm_boot2,[],length(fns),length(filled_ind));
fnind = 0;
for fn= fns
    fnind = fnind+1;
        %predictions
         if length(whichside)>1
            Preds(:,:,1) = Xm_boot(:,fnind,1:length(filled_ind));
            Preds(:,:,2) = Xm_boot2(:,fnind,1:length(filled_ind));
        elseif whichside ==1
            Preds(:,:,whichside) = Xm_boot(:,fnind,1:length(filled_ind));
        else 
            Preds(:,:,whichside) = Xm_boot2(:,fnind,1:length(filled_ind));
        end 
        
        scalePreds(fnind,:,:,whichside)=reshape(Preds(:,:,whichside),size(scalePreds(fnind,:,:,whichside)));
        
        for index = 1:length(filled_ind)
            if length(whichside)>1
                %Truth: lower half (preds1)
                
                Truthscale(:,index,1) = Xscale(row(index),col(index),:,Tempind);
                Truthscale(:,index,2) = Xscale(col(index),row(index),:,Tempind);

            elseif whichside==1
                %Truth: lower half (preds1)
                tempX = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
                filled_ind = find(~isnan(tempX));
                [row,col]= find(~isnan(tempX));
                Truthscale(:,index,1) = Xscale(row(index),col(index),:,Tempind);
            elseif whichside==2
                %top half truth
                
                Truthscale(:,index,2) = Xscale(col(index),row(index),:,Tempind);
            end 
        end 
     for c =1:size(X,3) 
        errors(fnind,c,:,whichside) = reshape(reshape(Truthscale(c,:,whichside),size(Preds(c,:,whichside)))-Preds(c,:,whichside), size(errors(fnind,c,:,whichside)));
        ARD(fnind,c,:,whichside) = reshape(abs(reshape(Truthscale(c,:,whichside),size(Preds(c,:,whichside)))-Preds(c,:,whichside))./reshape(Truthscale(c,:,whichside),size(Preds(c,:,whichside))),size(errors(fnind,c,:,whichside)));
        mse(fnind,c) = sum(sum(errors(fnind,c,:,whichside).^2))/length(errors(fnind,c,:,whichside))/2;
        wmse(fnind,c) = find_wmse_error(errors(fnind,c,:,whichside), length(filled_ind')*2);
        aard(fnind,c) = sum(sum(abs(ARD(fnind,c,:,whichside))))/length(ARD(fnind,c,:,whichside))/2;
        perclessscale(fnind,c,1) = sum(sum((abs(ARD(fnind,c,:,:))<0.05)))/length(filled_ind)/2*100;
        perclessscale(fnind,c,2) = sum(sum((abs(ARD(fnind,c,:,:))<0.15)))/length(filled_ind)/2*100;
        perclessscale(fnind,c,3) = sum(sum((abs(ARD(fnind,c,:,:))<0.25)))/length(filled_ind)/2*100;
        perclessscale(fnind,c,4) = sum(sum((abs(ARD(fnind,c,:,:))<0.50)))/length(filled_ind)/2*100;
        perclessscale(fnind,c,5) = sum(sum((abs(ARD(fnind,c,:,:))<0.75)))/length(filled_ind)/2*100;
    end  
end 

%%  He preds
load(filename)
%
mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 
Temperature = Temps(Tind);
fnind=0;
load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
[~,mixunind] = find(ismember(mixture',mixtures(filled_ind,:), 'rows')); % hopefully this extracts what we want 
heunifacpred = he(5:5:96, mixunind);
mixtureuni = mixture;
clear errorUNI
clear errorHE aardHE
for fnind =1:length(fns)

    rSign = 1;
    rScale = fnind;
    
    
    %hepred = signPreds(rSign,:,:,:).*exp(signPreds(rSign,:,:,:).*scalePreds(rScale,:,:,:));
    
    hepred = Truthsign(:,:,whichside).*exp(reshape(scalePreds(rScale,:,:,whichside),size(Truthsign(:,:,whichside))).*Truthsign(:,:,whichside));
    errorUNI(:,:,whichside) = hepred - heunifacpred; 
    errorUNIsystem = mean(errorUNI,1);
    %hepred = reshape(hepred,size(X,3),[],length(whichside));
    %hepred(:,:,2) = flip(hepred(:,:,2));
    %postprocess the predictions 
    if postprocess ==1
        [heprednew,errorfit] = postprocesshe(hepred,concentrations);
        hepred(:,:,1) = heprednew;
    end 
    
    Truth = zeros(size(hepred));
    Xtemp = tril(Xs(:,:,1,Tempind),-1)+triu(nan(size(Xs(:,:,1,Tempind))));
    Xtemp(Xtemp==0)=nan;
   [row,col,elements]=find(~isnan(Xtemp));
   
    for counter = 1:size(row,1)
        
        if length(whichside)>1
            Truth(:,counter,1) = X(row(counter),col(counter),:,Tempind);
            Truth(:,counter,2) = X(col(counter),row(counter),:,Tempind);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        elseif whichside ==1
            Truth(:,counter,1) = X(row(counter),col(counter),:,Tempind);
        elseif whichside ==2
            Truth(:,counter,2) = X(col(counter),row(counter),:,Tempind);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        end 
        
    end 
    Truth(isnan(Truth)) =0;
    indicespred = find(~isnan(hepred(1,:,1)));
    if whichside ==2
        errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred));
    else 
        errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred,whichside));
    end 
    errorHE(isnan(errorHE)) = 0;
    errorHEsystem = mean(errorHE,1);
    hepred = hepred(:,indicespred,:);
    aardHE(:,:,whichside) = abs(errorHE(:,:,whichside))./abs(reshape(Truth(:,indicespred,whichside),size(errorHE(:,:,whichside))));
    %aardHE(:,:,2) = abs(errorHE(:,:,2))./abs(reshape(Truth(:,indicespred,2),size(errorHE(:,:,2))));

    for c = 1:size(X,3)
        smsec(fnind,c) = sqrt(sum(sum(errorHE(c,:,whichside).^2))/length(indicespred));
        wsmsec(fnind,c) = sqrt(find_wmse_error(errorHE(c,:,whichside),length(indicespred)*2));
        aardc(fnind,c) = sum(sum(abs(aardHE(c,:,whichside))))/length(aardHE(c,:,1));
        percless(fnind,c,1) = sum(sum((abs(aardHE(c,:,whichside))<0.05)))/length(indicespred)/2*100;
        percless(fnind,c,2) = sum(sum((abs(aardHE(c,:,whichside))<0.15)))/length(indicespred)/2*100;
        percless(fnind,c,3) = sum(sum((abs(aardHE(c,:,whichside))<0.25)))/length(indicespred)/2*100;
        percless(fnind,c,4) = sum(sum((abs(aardHE(c,:,whichside))<0.50)))/length(indicespred)/2*100;
        percless(fnind,c,5) = sum(sum((abs(aardHE(c,:,whichside))<0.75)))/length(indicespred)/2*100;
    end 
    wsmse(fnind) = sqrt(find_wmse_error(errorHE,length(indicespred)));
    smse(fnind) = sqrt(sum(sum(sum(errorHE.^2)))/prod(size(errorHE)));
end
wsmse

%% Pure he preds
clc
clear
load('3wayHOSVD-threshold60-fn3=4-Small-concslices-LOOCV-Xnone-maxiter=20000-Temp=318.15-fillmethod=avr-20-Dec-2022.mat')
load(filename)
Tempind =4;
Temperature = Temps(Tempind); 
whichside = 1:2;
concen = 1;
postprocess = 0;
concentrations = conc_interval;
%load('3wayHOSVD-threshold60-centre+scale-Small-concslices-LOOCV-Xnone-maxiter=20000-Temp=318.15-fillmethod=avr-09-Dec-2022.mat')
tempX = tril(Xs(:,:,1,Tempind),-1)+triu(nan(size(Xs(:,:,1,Tempind))));
[row,col] = find(~isnan(tempX));
filled_ind =  find(~isnan(tempX));
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end

load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
[mixunind,~] = (ismember(mixture',mixtures(filled_ind,:), 'rows')); % hopefully this extracts what we want 
heunifacpred = he(5:5:96, mixunind);
mixtureuni = mixture;
clear errorUNI
clear errorHE aardHE
% need to find which fns were not used 
for fnind =1:length(fns)
    hepred = [Xm_boot(fnind,:,:); Xm_boot2(fnind,:,:)]; 
    hepred = permute(hepred, [3 2 1]);
    Truth = zeros(size(hepred));
    Xtemp = tril(Xs(:,:,1,Tempind),-1)+triu(nan(size(Xs(:,:,1,Tempind))));
    Xtemp(Xtemp==0)=nan;
    rSign = 1;
    rScale = fnind;
    
 
    if postprocess ==1
        [heprednew,errorfit] = postprocesshe(hepred(:,:,whichside),concentrations);
        hepred(:,:,whichside) = heprednew;
    end 
    
    Truth = zeros(size(hepred));
    Xtemp = tril(Xs(:,:,1,Tempind),-1)+triu(nan(size(Xs(:,:,1,Tempind))));
    Xtemp(Xtemp==0)=nan;
   [row,col,elements]=find(~isnan(Xtemp));
   
    for counter = 1:size(row,1)
        
        if length(whichside)>1
            Truth(:,counter,1) = X(row(counter),col(counter),:,Tempind);
            Truth(:,counter,2) = X(col(counter),row(counter),:,Tempind);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        elseif whichside ==1
            Truth(:,counter,1) = X(row(counter),col(counter),:,Tempind);
        elseif whichside ==2
            Truth(:,counter,2) = X(col(counter),row(counter),:,Tempind);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        end 
        
    end 
    Truth(isnan(Truth)) =0;
    indicespred = find(~isnan(hepred(1,:,1)));
    
    errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred,whichside));
    errorHE(errorHE==0)=nan; 
    errorUNI(:,:) = Truth(:,:,1) - heunifacpred; 
    errorUNIsystem = mean(errorUNI,1);
    errorHE(isnan(errorHE)) = 0;
    errorHEsystem = mean(errorHE,1);
    hepred = hepred(:,indicespred,:);
    aardHE(:,:,whichside) = abs(errorHE(:,:,whichside))./abs(reshape(Truth(:,indicespred,whichside),size(errorHE(:,:,whichside))));
    %aardHE(:,:,2) = abs(errorHE(:,:,2))./abs(reshape(Truth(:,indicespred,2),size(errorHE(:,:,2))));

    for c = 1:size(X,3)
        smsec(fnind,c) = sqrt(sum(sum(errorHE(c,:,whichside).^2))/length(indicespred));
        wsmsec(fnind,c) = sqrt(find_wmse_error(errorHE(c,:,whichside),length(indicespred)*2));
        aardc(fnind,c) = sum(sum(abs(aardHE(c,:,whichside))))/length(aardHE(c,:,1));
        percless(fnind,c,1) = sum(sum((abs(aardHE(c,:,whichside))<0.05)))/length(indicespred)/2*100;
        percless(fnind,c,2) = sum(sum((abs(aardHE(c,:,whichside))<0.15)))/length(indicespred)/2*100;
        percless(fnind,c,3) = sum(sum((abs(aardHE(c,:,whichside))<0.25)))/length(indicespred)/2*100;
        percless(fnind,c,4) = sum(sum((abs(aardHE(c,:,whichside))<0.50)))/length(indicespred)/2*100;
        percless(fnind,c,5) = sum(sum((abs(aardHE(c,:,whichside))<0.75)))/length(indicespred)/2*100;
    end 
    wsmse(fnind) = sqrt(find_wmse_error(errorHE,length(indicespred)));
    smse(fnind) = sqrt(sum(sum(sum(errorHE.^2)))/prod(size(errorHE)));
end
wsmse

%%
f1=figure(1);
f1.Position = [10,10,400, 300];
load('colorblind_colormap.mat')
clf
index1plot=1;
semilogy(fns(index1plot:end),smse(index1plot:end),'.','Color', colorblind(1,:), 'MarkerSize',14)
hold on 
semilogy(fns(index1plot:end),wsmse(index1plot:end),'o','Color', colorblind(2,:), 'MarkerSize',6)

ylabel('Error (J/mol)')
xlabel('Number of factors')
legend('SMSE','wSMSE')
%%
clf
index1plot=1;
%choose mse or wmse
plot(fns(index1plot:end),wsmse(index1plot:end),'.','Color', colorblind(1,:),'MarkerSize',14)
ylabel('wSMSE (J/mol)')

xlabel('Number of factors')
%%
plot(fns(index1plot:end),wsmse(index1plot:end),'o','Color', colorblind(2,:),'MarkerSize',6)
ylabel('wSMSE (J/mol)')

xlabel('Number of factors')

%% Define mixture you wish to predict
figure(2)
clf

p = 0;

load('colorblind_colormap.mat')
%info we might want later 
func_group = {1,2,3:5, [6,7,11], 8:9, 10, 12:17,18:20, 21, 22};
label_func = {'Alkane', 'Primary alcohol', 'Other alcohol','Cycloalkanes', 'Ketone', 'Alkene', 'Ester', 'Amine','Acid','Aldehyde'};

%functional groups

%load experimental data from the array used in the data - should work if the array has been extended as well   
load('HE4wayArrayPolySmall4', 'conc_original', 'HE_original', 'mixtureT') 
mixture_exp = mixtureT;

concentrations = 0.05:0.05:0.95;

for counter =1:16
    
    index = counter+p*16;
    mixpredind = filled_ind(index);%change the mixture number here
    mixpred = mixtures(mixpredind,:); % might need to change this depending on the results 
    %match mixpred to mixture_exp 
    [~,indexpred] = ismember(mixpred,mixture_exp,'rows');
    
    conc_exp = conc_original(:,indexpred, Tempind);
    heexp = HE_original(:,indexpred, Tempind);
    
    disp(mixpred)
    
     subplot(4,4,counter)
    plot(conc_exp*100, heexp,'^','Color', colorblind(6,:), 'MarkerSize',4,'MarkerFaceColor',colorblind(6,:))
    hold on
    
    %predictions made by the model 
    if length(whichside)>1
        heplot1 = (hepred(:,counter+p*16,1));
        heplot2 = (hepred(:,counter+p*16,2));
        plot(concentrations*100, flip(heplot1),'.','Color', colorblind(8,:), 'MarkerSize',12)      
        hold on
        if postprocess == 0
            plot(concentrations*100, (heplot2) ,'.', 'Color',colorblind(1,:), 'MarkerSize',12)
            hold on
        end 
    else
        heplot1 = (hepred(:,counter+p*16,whichside));
        plot(concentrations*100, heplot1,'.', 'Color',colorblind(8,:), 'MarkerSize',12)
        hold on 
    end 
    
    %extract experimental predictions 
    
    
    %extract unifac predictions 
    load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
    [~,mixunind] = ismember(mixpred,mixture', 'rows');
    if mixunind ~=0
        heuni = he(:,mixunind);
        plot(conc_interval*100, (heuni),'k--')
    end 
    
    title(num2str(mixpred))% work on making this words - dictionary of compounds and mixture numeric? 
    %title(strcat(label_func(find(strcmp(func_group, {(mixpred(1,1))}))), ' ', num2str(mixpred(1,2)), '+', label_func(find(strcmp(func_group, {num2str(mixpred(1,3))}))), ' ',num2str(mixpred(1,4))))) 
    
    hold off
end 
xlabel('Composition of compound 1 (%)')
ylabel('Excess enthalpy (J/mol)')
if length(whichside)>1 && postprocess ==0
    legend('Experimental','Predictions left', 'Predictions right','UNIFAC (Do)')
else
    legend('Experimental','Predictions',  'UNIFAC (Do)')
end 
%% Errors per type of mixtures 

%Per type
dim = size(Truth);
type = zeros(1,dim(2));
for counter = 1:dim(2) % for all systems 
    %extract the data needed 
    temp = reshape(Truth(:,counter,whichside),dim(1),[]);
    maxVal = max(temp);
    minVal = min(temp);
    if length(whichside)>1
        maxVal = maxVal(1);
        minVal = minVal(1);
    end 
    if maxVal<0 || minVal<0 
        % passes through the origin or is negative 
        if sign(maxVal)~=sign(minVal) %passes through origin, only one is negative 
            type(counter) = 1;
        else 
            type(counter) = 2;
        end 
    else 
        %positive curve 
        if maxVal > 1000 % adjustable ceilings for this classification
            type(counter) = 3; % very positive 
        elseif maxVal >200
            type(counter) = 4; % moderately positive 
        else 
            type(counter) = 5; % less than 200 at the maximum 
        end 
    end 
end
%errors per type 

for counter = 1:5
    temp = errorHE(:,find(type==counter),:);
    temp2 = aardHE(:,find(type==counter),:);
    errortype{counter} =errorHE(:,find(type==counter),:);
    smsetype(counter) = sqrt(sum(sum(sum(temp.^2)))/prod(size(temp)));
    wsmsetype(counter) = sqrt(find_wmse_error(temp,prod(size(temp))));
    aardtype(counter) = sum(sum(sum(abs(temp2))))/prod(size(temp2));
    nomix(counter) = length(find(type==counter));
end 
tabletypes = [1:5; nomix;smsetype; wsmsetype; aardtype];

%% Errors per type of functional groups in the mixture 
figure(3)
clf
func_group = {1,2,3:5, [6,7,11], 8:9, 10, 12:17,18:20, 21, 22};
label_func = {'Alkane', 'Primary alcohol', 'Other alcohol','Cycloalkanes', 'Ketone', 'Alkene', 'Ester', 'Amine','Acid','Aldehyde'};

err_funcgroup = cell(length(func_group),length(func_group));
ard_funcgroup = cell(length(func_group),length(func_group));
smse_funcgroup = zeros(length(func_group),length(func_group));
wsmse_funcgroup = zeros(length(func_group),length(func_group));
AARD_funcgroup = zeros(length(func_group),length(func_group));
no_mix = zeros(length(func_group),length(func_group));

%extracting the correct mixtures, if necessary 

mixtureT = intersect(mixtures,mixtureT,'rows');
mixtureT = mixtureT(1:length(filled_ind),:);
clf
for func = 1:length(func_group)
    for func2 = 1:length(func_group)
        %extract arrays of groups 
        funcgroup1 = func_group{func};
        funcgroup2 = func_group{func};
        %find mixtures with this combination
        index = find(mixtureT(:,1)==func);
        index2 = find(mixtureT(:,3)==func2);
        indices = intersect(index,index2);
        %export errors and AARDs
        err_funcgroup{func,func2} = errorHE(:,indices,:); 
        errtemp = errorHE(:,indices,:);
        smse_funcgroup(func,func2) = sqrt(sum(errtemp(:).^2)/length(errtemp(:)));
        wsmse_funcgroup(func,func2) = find_wmse_error(errtemp(:),length(errtemp(:)));
        ard_funcgroup{func,func2} = aardHE(:,indices,:);
        ardtemp = aardHE(:,indices,:);
        AARD_funcgroup(func,func2) = sum(ardtemp(:))/length(ardtemp(:));  
        no_mix(func,func2) = length(ardtemp(:)/9);
    end 
end 

% Heatmap of the errors per functional group 

SMSEplot = smse_funcgroup;
AARDplot=AARD_funcgroup.*100; 
%create a table 
[row,col] = find(~isnan(AARDplot));
for counter = 1:length(row)
    funcgroup1tbl{counter} = label_func{row(counter)};
    funcgroup2tbl{counter} = label_func{col(counter)};
    smse_tbl(counter) = SMSEplot(row(counter), col(counter));
    aard_tbl(counter) = AARDplot(row(counter),col(counter));
    wsmsetbl(counter) = wsmse_funcgroup(row(counter),col(counter));
end 

% heatmap 
tbl = table(funcgroup1tbl', funcgroup2tbl',smse_tbl',aard_tbl',wsmsetbl', 'VariableNames', ["Functional group 1", "Functional group 2", "SMSE (J/mol)","AARD (%)","wSMSE (J/mol)" ]);
% heatmap 
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'SMSE (J/mol)', 'ColorMethod', 'none');
h.Colormap=winter;
h.FontSize = 14;

%% Parity plot and error distribution plots 
%Import all the UNIFAC predictions 
load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
for i = 1:length(filled_ind)
    mixpred = mixtureT(i,:);
    [~,mixunind] = ismember(mixpred,mixture', 'rows');
    if mixunind ~=0
        heuniplot(:,i) = he(5:5:96,mixunind); % check these indices 
        indexuni(i)=1;
    else 
        indexuni(i)=0;
    end 
end 
f4 = figure(4);
f4.Position=[10 10 400 300];
clf 
load('colorblind_colormap.mat')

%choose between plot and loglog 
plotx =Truth(:,indicespred,whichside);
ploty = hepred(:,indicespred,whichside);
plot(plotx(:),ploty(:),'.','Color', colorblind(1,:), 'MarkerSize',8)
hold on 
plot(plotx(:,logical(indexuni),1), heuniplot,'.','Color', 'k', 'MarkerSize',8)
%prediction interval
PIx = min(plotx(:)):100:max(plotx(:));
syx = std(plotx(:));
hi = 1/numel(plotx(:)) + (PIx-mean(PIx)).^2/sum((PIx-mean(PIx)).^2);
PIy1 = PIx - 1.96.*syx.*sqrt(1+hi);
PIy2 = PIx + 1.96.*syx.*sqrt(1+hi);
plot(PIx,PIy1, '--', 'Color', colorblind(7,:), 'LineWidth', 1) %make this pink 
plot(PIx,PIy2, '--', 'Color', colorblind(7,:), 'LineWidth', 1) %make this pink 
%y=x line 
plot([min(plotx(:)),max(plotx(:))], [min(plotx(:)),max(plotx(:))],'-', 'Color',colorblind(4,:), 'LineWidth', 1.2)
 

hold off 
xlabel('Experimental (J/mol)','FontSize',13, 'Color', 'k')
ylabel('Prediction (J/mol)','FontSize',13, 'Color', 'k')
if plotuni
    legend('MCM', 'UNIFAC','PI (95%)', 'FontSize',8,'Location', 'northwest');
else 
    legend('MCM','PI (95%)', 'FontSize',8,'Location', 'northwest');
end 

%error distribution 
f5 = figure(5);
f5.Position = [10 20 400 300];
clf
%choose the number of bins 
%ploterrhist(errorHE(:), 'bins', 20)
histogram(errorHE(:), 20, 'FaceColor', colorblind(6,:), 'BinLimits', [-1000 1000])
ylabel('Instances','FontSize',13, 'Color', 'k')
xlabel('MCM error (J/mol)', 'FontSize',13, 'Color', 'k')

%% 3way histogram 
f6 = figure(6);
f6.Position = [10 20 600 400];
clf
edges = {-1000:100:1000; 0:5:100}; % bin edges
data = [errorHE(:) abs(aardHE(:)).*100]; 
hist3(data, 'Edges', edges,'CdataMode','auto','FaceColor','interp')
zlabel('Instances','FontSize',13, 'Color', 'k')
xlabel('MCM error (J/mol)', 'FontSize',13, 'Color', 'k')
ylabel('MCM AARD (%)' ,'FontSize', 13, 'Color', 'k')
colorbar
colormap('parula')

%% Error plot within the 3-way arrays 
f7 = figure(7);
f7.Position = [10 20 800 500];
clf
TempIndex = 1;
errorplot = nan(size(X(:,:,:,TempIndex)));
for i =1:length(col)
    errorplot(row(i),col(i),:) = errorHE(:,i,1);
    errorplot(col(i),row(i),:) = errorHE(:,i,2);
end 
[X,Y,Z] = meshgrid(1:size(errorplot,1), 1:size(errorplot,1), 5:5:95);
xslice = 1:1:size(errorplot,1);    % location of y-z planes
yslice = 1:1:size(errorplot,1);     % location of x-z plane
zslice = 5:5:95;         % location of x-y planes
slice(X,Y,Z,errorplot,xslice,yslice,zslice)
c = colorbar;
c.Label.String = 'Error (J/mol)';
c.Label.FontSize = 12;
%% Classification of correct predictions per system 
% colours of data points will show which one is better 
%started 04/11/2022 - check if it works 
%scatter(x,y,pointsize,color)

leftcorrect = errorHEsystem(1,:,1)<errorHEsystem(1,:,2);% logical output 
errorHEbest = errorHEsystem(1,:,1).*leftcorrect + errorHEsystem(1,:,2).*(1-leftcorrect);
unifaccorrect = errorHEbest>(errorUNIsystem(1,:,1).*leftcorrect+ errorUNIsystem(1,:,1).*(1-leftcorrect)); %logical output 
numleftcorrect = sum(find(leftcorrect));
numunicorrect = sum(find(unifaccorrect));
%figure out x and y later 
scatter(1:length(filled_ind),1:length(filled_ind),[],leftcorrect, 'filled')
scatter(1:length(filled_ind),1:length(filled_ind),[],unifaccorrect, 'filled')

%Classification by machine learning method 
%% Classification of correct predictions per functional group 
%fix all of this before using - 04/11/2022 - does not do anything 

for func = 1:length(func_group)
    for func2 = 1:length(func_group)
        %extract arrays of groups 
        funcgroup1 = func_group{func};
        funcgroup2 = func_group{func};
        %find mixtures with this combination
        index = find(mixtureT(:,1)==func);
        index2 = find(mixtureT(:,3)==func2);
        indices = intersect(index,index2);
        %export errors and AARDs
        left_funcgroup{func,func2} = sum(leftcorrect(indices)); 
        unifac_funcgroup{func,func2} = sum(unifaccorrect(indices));
        ardtemp = aardHE(:,indices,:); 
        no_mix(func,func2) = length(ardtemp(:)/19);
    end 
end 
 
%create a table

[row,col] = find(~isnan(AARDplot));
% Find which side was more correct and whether UNIFAC or the model made better predictions 
for counter = 1:length(row)
    funcgroup1tbl{counter} = label_func{row(counter)};
    funcgroup2tbl{counter} = label_func{col(counter)};
    leftcorrect_tbl(counter) = left_funcgroup{row(counter),col(counter)};%left side correct = 1
    unifaccorrect_tbl(counter) =unifac_funcgroup{row(counter),col(counter)};
    
end 

% heatmap 
tbl = table(funcgroup1tbl', funcgroup2tbl',leftcorrect_tbl',unifaccorrect_tbl', 'VariableNames', ["Functional group 1", "Functional group 2", "Left side correct ","UNIFAC correct" ]);
% heatmap - color = % correct per functional group 
figure(5)
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'UNIFAC correct', 'ColorMethod', 'none');
h.Colormap=winter;
h.FontSize = 14;
figure(6)
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'Left side correct', 'ColorMethod', 'none');
h.Colormap=winter;
h.FontSize = 14;


%% functions 
function [heprednew,uncertainty] = postprocesshe(hepred,concentrations)
    heprednew = zeros(size(hepred,1),size(hepred,2));
    
    uncertainty = zeros(length(size(hepred,2)),length(concentrations));
    threshold = 0.5*(size(hepred,1)*2-1)/sqrt(size(hepred,1)*2);

    for index = 1:(size(hepred,2))
        temppred = hepred(:,index,:);
        xvar = [concentrations,concentrations];
        temppred = temppred(:);
        %remove abnormal values and zeros 
        temppred(abs(temppred)>1e12)=nan;
        temppred(temppred==0)=nan;
        temppred = temppred(~isnan(temppred));
        xvar = xvar(~isnan(temppred));
        % Use the threshold to remove outliers 
        tempscore = (temppred-mean(temppred))/std(temppred);
        temppred(abs(tempscore)>threshold)=nan;
        temppred = temppred(~isnan(temppred));
        xvar = xvar(~isnan(temppred));
        xvar = [xvar,0,1];
        temppred = [temppred; 0;0];
        %disp(xvar)
        %disp(temppred)
        %fit a polynomial 
        [p,S,mu] = polyfit(xvar, temppred,3);
        % create polynomial predictions 
        [temp,uncertainty(index,:)] = polyval(p,concentrations,S,mu);
        heprednew(:,index) = temp;
    end 
    
end 
function [X_removed]=remove_nan2(X)
    %check for rows and columns containing only nan values in a matrix and remove them 
    % Input 
    % X = matrix with rows or columns containing nan values
    % Output 
    % X_removed = matrix without rows or columns containing only nan values
    % 
    row_nan = find(all(isnan(X),2));
    col_nan = find(all(isnan(X),1));
    X_removed = X;
    if row_nan
        X_removed(row_nan,:)=[];
    end 
    if col_nan 
        X_removed(:,col_nan)=[];
    end
end 