%% Francesca Middleton 27/09/2022

%discovering the percent correct prediction of the sign predictions 

clc
clear
load('colorblind_colormap.mat')
Temperature=298.15;
concen = 1;
ways = 2;
%if only one side of the 2-way array is considered or both (1:2)
whichside =1:2;
postprocess = 0;
load('2waySVD-33comps-threshold60par-LOOCV-Xnone-T=293.15-fillmethod=avr-19-Dec-2022.mat')
%load("2waySVD-71comps-threshold60par-LOOCV-Xscale-T=298.15-fillmethod=reg-19-Nov-2022.mat")
%load("2waySVD-27comps-addingfrom20-threshold60par-LOOCV-Xscale-T=298.15-fillmethod=reg-22-Nov-2022.mat")
%load('2waySVD-27comps-noadding-threshold60par-LOOCV-Xscale-T=298.15-fillmethod=reg-16-Nov-2022.mat')
%% Xsign
load("2waySVD-27comps-addingfrom27-threshold60par-LOOCV-Xsign-T=298.15-fillmethod=reg-23-Nov-2022.mat")

Xtemp = Xs;
Xtemp(X==0)=nan;
tempX = Xtemp(:,:,1);
indend = 00; %this is valid for when compounds were added
tempX(1:indend,1:indend) = nan;
tempX = tril(tempX,-1)+triu(nan(size(tempX)));
[row,col] = find(~isnan(tempX));
filled_ind = find(~isnan(tempX));
signPrediction = sign(Xm_boot);
signPrediction(signPrediction ==0)=nan;
indicesSign = find(~isnan(signPrediction(1,1,:)));
signPrediction = signPrediction(:,:,indicesSign);
errors = zeros(length(fns),size(X,3),length(filled_ind),2);
signPreds = zeros(length(fns),size(X,3),length(filled_ind),2);
correct = zeros(length(fns),size(X,3));
%find errors and how many were correctly predicted 

Preds = [Xm_boot(1,1,1:length(filled_ind))];%; Xm_boot2(1,1,:)];
%Preds = reshape(Preds,2,[]);
indiceskeep = find(Preds(1,:)); % not equal to zero 
fnind = 0;
for fn = 1:length(fns)
    fnind = fnind+1;
    for c = 1:size(X,3)

        if length(whichside)>1
            Preds(1,:) = Xm_boot(c,fn,1:length(filled_ind));
            Preds(2,:) = Xm_boot2(c,fn,1:length(filled_ind));
        elseif whichside ==1
            Preds(1,:) = Xm_boot(c,fn,1:length(filled_ind));
        else 
            Preds(1,:) = Xm_boot2(c,fn,1:length(filled_ind));
        end 
        
        Preds = sign(Preds(:, indiceskeep));
        signPreds(fn,c,:,whichside) =reshape(Preds,size(signPreds(1,1,:,whichside))); 
        
        if length(whichside)>1
            %Truth: lower half (preds1)
            
            Truthtemp = Xs(:,:,c);

            Truthsign(fnind,c,:,1) = Truthtemp(filled_ind);
            %top half truth
            tempX = (triu(Xs(:,:,c),1)+tril(nan(size(Xs(:,:,c)))))';
            filled_ind = find(~isnan(tempX));
            Truthsign(fnind,c,:,2) = (tempX(filled_ind));
        elseif whichside==1
            %Truth: lower half (preds1)
            tempX = tril(Xs(:,:,c),-1)+triu(nan(size(Xs(:,:,c))));
            filled_ind = find(~isnan(tempX));
            Truthtemp = Xs(:,:,c);

            Truthsign(fnind,c,:,1) = Truthtemp(filled_ind);
        elseif whichside==2
            %top half truth
            tempX = (triu(Xs(:,:,c),1)+tril(nan(size(Xs(:,:,c)))))';
            filled_ind = find(~isnan(tempX));
            Truthsign(fnind,c,:,2) = (tempX(filled_ind));
        end 
        
        errors(fn,c,:,whichside) =reshape(reshape(Truthsign(fnind,c,:,whichside),size(Preds))-Preds, size(errors(1,1,:,whichside)));
        correct(fn,c)= sum(sum(errors(fn,c,:,whichside)==0))/length(Preds)/(length(whichside))*100;
        
    end 
end 

%% If no sign predictions exist

%finds the true sign values
X(X==0)=nan;
for counter =1:size(X,3)
    Xnew(:,:,counter) = remove_nan2(X(:,:,counter));
end 
dim = size(Xs);
X = Xnew;
Xsign = sign(X);
concentrations = conc_interval;

for j = 1:dim(3)
    Xtemp = Xs(:,:,j);
    Xtemp(1:1+size(Xtemp,1):end) = nan; % diagonals are nan
    Xs(:,:,j) = Xtemp;
end 
tempX = Xs(:,:,1);

indend = 20; %this is valid for when compounds were added
tempX(1:indend,1:indend) = nan;
tempX = tril(tempX,-1)+triu(nan(size(tempX)));
[row,col] = find(~isnan(tempX));
filled_ind = find(~isnan(tempX));
         

for index = 1:length(filled_ind)
    Truthsign(1,:,index,1) = Xsign(row(index),col(index),:);
    Truthsign(1,:,index,2) = (Xsign(col(index),row(index),:)); 
end 
load(filename)
mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 


%% Xscale
%load('2waySVD-27comps-noadding-threshold60par-LOOCV-Xscale-T=298.15-fillmethod=reg-16-Nov-2022.mat', 'Xm_boot', 'Xs', 'Xm_boot2','fns', 'filename')
load(filename)
mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 

concentrations = 0.05:0.05:0.95;
clear Preds
scalePrediction = Xm_boot;
scalePrediction(scalePrediction ==0)=nan;
indicesScale = find(~isnan(scalePrediction(1,1,:)));

Xtemp = Xs;
Xtemp(X==0)=nan;
tempX = Xtemp(:,:,1);
tempX(1:indend,1:indend) = nan;
tempX = tril(tempX,-1)+triu(nan(size(tempX)));
[row,col] = find(~isnan(tempX));
filled_ind = find(~isnan(tempX));

X(X==0)=nan;
for counter =1:size(X,3)
    Xnew(:,:,counter) = remove_nan2(X(:,:,counter));
end 

X = Xnew;
Xscale = log(sign(X).*(X)).*sign(X);
scalePreds = zeros(length(fns),size(X,3),length(filled_ind),2);
mse = zeros(length(fns),size(X,3));
wmse = zeros(length(fns),size(X,3));
aard = zeros(length(fns),size(X,3));
ARD = zeros(length(fns),size(X,3),length(filled_ind),2);
errors = zeros(length(fns),size(X,3),length(filled_ind),2);

Preds = zeros(length(concentrations),length(filled_ind),length(whichside));
%Preds = reshape(Preds,2,[]);

for fn = 1:length(fns)
    
        %predictions
         if length(whichside)>1
            Preds(:,:,1) = Xm_boot(:,fn,1:length(filled_ind));
            Preds(:,:,2) = Xm_boot2(:,fn,1:length(filled_ind));
        elseif whichside ==1
            Preds(:,:,whichside) = Xm_boot(:,fn,1:length(filled_ind));
        else 
            Preds(:,:,whichside) = Xm_boot2(:,fn,1:length(filled_ind));
        end 
        
        scalePreds(fn,:,:,whichside)=reshape(Preds(:,:,whichside),size(scalePreds(fn,:,:,whichside)));
        
        for index = 1:length(filled_ind)
            if length(whichside)>1
                %Truth: lower half (preds1)
                
                Truthscale(:,index,1) = Xscale(row(index),col(index),:);
                Truthscale(:,index,2) = Xscale(col(index),row(index),:);

            elseif whichside==1
                %Truth: lower half (preds1)
                tempX = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
                filled_ind = find(~isnan(tempX));
                [row,col]= find(~isnan(tempX));
                Truthscale(:,index,1) = Xscale(row(index),col(index),:);
            elseif whichside==2
                %top half truth
                tempX = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
                filled_ind = find(~isnan(tempX));
                [row,col]= find(~isnan(tempX));
                Truthscale(:,index,2) = Xscale(col(index),row(index),:);
            end 
        end 
     for c =1:size(X,3)   
        errors(fn,c,:,whichside) = reshape(reshape(Truthscale(c,:,whichside),size(Preds(c,:,whichside)))-Preds(c,:,whichside), size(errors(fn,c,:,whichside)));
        ARD(fn,c,:,whichside) = reshape(abs(reshape(Truthscale(c,:,whichside),size(Preds(c,:,whichside)))-Preds(c,:,whichside))./reshape(Truthscale(c,:,whichside),size(Preds(c,:,whichside))),size(errors(fn,c,:,whichside)));
        mse(fn,c) = sum(sum(errors(fn,c,:,whichside).^2))/length(errors(fn,c,:,whichside))/2;
        wmse(fn,c) = find_wmse_error(errors(fn,c,:,whichside), length(filled_ind')*2);
        aard(fn,c) = sum(sum(abs(ARD(fn,c,:,whichside))))/length(ARD(fn,c,:,whichside))/2;
        perclessscale(fn,c,1) = sum(sum((abs(ARD(fn,c,:,:))<0.05)))/length(filled_ind)/2*100;
        perclessscale(fn,c,2) = sum(sum((abs(ARD(fn,c,:,:))<0.15)))/length(filled_ind)/2*100;
        perclessscale(fn,c,3) = sum(sum((abs(ARD(fn,c,:,:))<0.25)))/length(filled_ind)/2*100;
        perclessscale(fn,c,4) = sum(sum((abs(ARD(fn,c,:,:))<0.50)))/length(filled_ind)/2*100;
        perclessscale(fn,c,5) = sum(sum((abs(ARD(fn,c,:,:))<0.75)))/length(filled_ind)/2*100;
    end  
end 


%%
load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
mixunind = find(ismember(mixture',mixtures(filled_ind,:), 'rows')); % hopefully this extracts what we want 
mixturepred = mixtures(filled_ind,:);
mixtureuni = mixture(:,mixunind)';
[La,Locb]=ismember(mixtureuni,mixturepred, 'rows');
%mixture
%mixunind = mixunind(Locb,:);
heunifacpred = he(5:5:96, mixunind);
 rSign = 1;
for fnind =19%1:length(fns)

   
    rScale = fnind;
    
    
    %hepred = signPreds(rSign,:,:,:).*exp(signPreds(rSign,:,:,:).*scalePreds(rScale,:,:,:));
    
    hepred = Truthsign(rSign,:,:,whichside).*exp(reshape(scalePreds(rScale,:,:,whichside),size(Truthsign(rSign,:,:,whichside))).*Truthsign(rSign,:,:,whichside));
    hepred = reshape(hepred, size(Truthsign,2), size(Truthsign,3), size(Truthsign,4));
    %hepred = reshape(hepred,size(X,3),[],length(whichside));
    %hepred(:,:,2) = flip(hepred(:,:,2));
    %postprocess the predictions 
    if postprocess ==1
        [heprednew,errorfit] = postprocesshe(hepred,concentrations);
        hepred(:,:,1) = heprednew;
    end 
    
    Truth = zeros(size(hepred));
    Xtemp = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
    Xtemp(Xtemp==0)=nan;
    
   %[row,col,elements]=find(~isnan(Xtemp));
   
    for counter = 1:size(row,1)
        
        if length(whichside)>1
            Truth(:,counter,1) = X(row(counter),col(counter),:);
            Truth(:,counter,2) = X(col(counter),row(counter),:);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        elseif whichside ==1
            Truth(:,counter,1) = X(row(counter),col(counter),:);
        elseif whichside ==2
            Truth(:,counter,2) = X(col(counter),row(counter),:);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        end 
        
    end 
    Truth(isnan(Truth)) =0;
    indicespred = find(~isnan(hepred(concen,:,1)));
    if whichside ==2
        errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred));
    else 
        errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred,whichside));
    end 
    errorHE(isnan(errorHE)) = 0;
    errorHEsystem{fnind} = mean(errorHE,1);
    errorHEcomp{fnind} = mean(errorHE,2);
    errorUNI(:,:) = Truth(:,1:size(heunifacpred,2),1) - heunifacpred; 
    errorUNIsystem = mean(errorUNI,1);
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
indexplot2 = length(fns);
semilogy(fns(index1plot:indexplot2),smse(index1plot:indexplot2),'.','Color', colorblind(1,:), 'MarkerSize',14)
hold on 
semilogy(fns(index1plot:indexplot2),wsmse(index1plot:indexplot2),'o','Color', colorblind(2,:), 'MarkerSize',6)
hold off
ylabel('Error (J/mol)')
xlabel('Number of factors')
legend('SMSE','wSMSE')
%%
clf
index1plot=5;
%choose mse or wmse
plot(fns(index1plot:end),smse(index1plot:end),'.','Color', colorblind(1,:),'MarkerSize',14)
ylabel('SMSE (J/mol)')

xlabel('Number of factors')
%%
clf
indexplot=5;
plot(fns(index1plot:end),wsmse(index1plot:end),'o','Color', colorblind(2,:),'MarkerSize',6)
ylabel('wSMSE (J/mol)')

xlabel('Number of factors')
%% Pure he predictions 
clc 
clear
Temperature = 313.15; 
whichside = 1:2;
concen = 1;
postprocess = 0;
load('2waySVD-38comps-threshold60par-LOOCV-Xnone-T=313.15-fillmethod=avr-19-Dec-2022.mat')
concentrations = conc_interval;
%load("2waySVD-20comps-addingfrom20-threshold60par-LOOCV-Xnone-T=298.15-fillmethod=avr-24-Nov-2022.mat")
load(filename, 'comps')
    mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 
Xtemp = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
Xtemp(Xtemp==0)=nan;
filled_ind = find(~isnan(Xtemp));
[row,col] =  find(~isnan(Xtemp));

load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
mixunind = find(ismember(mixture',mixtures(filled_ind,:), 'rows')); % hopefully this extracts what we want 
mixturepred = mixtures(filled_ind,:);
mixunind2 = find(ismember(mixture', mixtures(filled_ind,[3,4,1,2]),'rows'));
mixuniInd = union(mixunind,mixunind2);
mixtureuni = mixture(:,mixuniInd)';
[La,Locb]=ismember(mixtureuni,mixturepred, 'rows');
%mixture
%mixunind = mixunind(Locb,:);
heunifacpred = he(5:5:96, mixuniInd);
heunifacpred(:, length(mixunind):length(mixuniInd)) = flip(heunifacpred(:, length(mixunind):length(mixuniInd)));
Xm_boot(Xm_boot==0)=nan;
Xm_boot = Xm_boot(~isnan(Xm_boot));
Xm_boot = reshape(Xm_boot, length(concentrations),length(fns),[]);
Xm_boot2(Xm_boot2==0)=nan;
Xm_boot2 = Xm_boot2(~isnan(Xm_boot2));
Xm_boot2 = reshape(Xm_boot2, length(concentrations),length(fns),[]);
[compoundnames, codesNames] = findnames([mixturepred(:,[1,2]); mixturepred(:,[3,4])]);

for fnind =1%1:length(fns)% 
    hepred = [Xm_boot(:,fnind,:) Xm_boot2(:,fnind,:)]; 
    hepred = permute(hepred, [1 3 2]);
    Truth = zeros(size(hepred));
    Xtemp = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
    Xtemp(Xtemp==0)=nan;

    %[row,col,elements]=find(~isnan(Xtemp));

    for counter = 1:size(row,1)

        if length(whichside)>1
            Truth(:,counter,1) = X(row(counter),col(counter),:);
            Truth(:,counter,2) = X(col(counter),row(counter),:);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        elseif whichside ==1
            Truth(:,counter,1) = X(row(counter),col(counter),:);
        elseif whichside ==2
            Truth(:,counter,2) = X(col(counter),row(counter),:);
            %Truth(:,i,2) = flip(Truth(:,i,2));
        end 

    end 
    Truth(isnan(Truth)) =0;
    indicespred = find(~isnan(hepred(concen,:,1)));
     
    errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred,whichside));
     
    errorHE(isnan(errorHE)) = 0;
    errorHEsystem{fnind} = mean(errorHE,1);
    errorHEcomp{fnind} = mean(errorHE,2);
    errorUNI(:,:) = Truth(:,1:size(heunifacpred,2),1) - heunifacpred; 
    errorUNIsystem = mean(errorUNI,1);
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

%% Define mixture you wish to predict
f2=figure(2);
f2.Position = [100,100,1000, 600];
load('colorblind_colormap.mat')
clf
p = 0;


%info we might want later 
func_group = {1,2,3:5, [6,7,11], 8:9, 10, 12:17,18:20, 21, 22};
label_func = {'Alkane', 'Primary alcohol', 'Other alcohol','Cycloalkanes', 'Ketone', 'Alkene', 'Ester', 'Amine','Acid','Aldehyde'};
load('colorblind_colormap.mat')
%functional groups
%filename might not be correct here
%load experimental data from the array used in the data - should work if the array has been extended as well   
load(filename, 'conc_original', 'HE_original', 'mixtureT') 
mixture_exp = mixtureT;

concentrations = 0.05:0.05:0.95;

for counter =1:16

    index = counter+p*16;
    mixpred = mixtureT(size(mixtureT,1)-length(filled_ind)+1*p+(counter),:);
    mixpredind = filled_ind(index);%change the mixture number here
    mixpred = mixtures(mixpredind,:); % might need to change this depending on the results 
    %match mixpred to mixture_exp 
    [~,indexpred] = ismember(mixpred,mixture_exp,'rows');
    
    if indexpred==0
        [~,indexpred] = ismember([mixpred(3:4) mixpred(1:2)], mixture_exp,'rows');
        if indexpred == 0
            indexpred = 1;
            disp('Mixture does not exist in experimental data')
            disp(mixpred)
        end 
    end 
    
    conc_exp = conc_original(:,indexpred);
    heexp = HE_original(:,indexpred);
    
    disp(mixpred)
    
    subplot(4,4,counter)
    plot(conc_exp*100, heexp,'^','Color', colorblind(6,:), 'MarkerSize',4,'MarkerFaceColor',colorblind(6,:))
    hold on
    
    %predictions made by the model 
    if length(whichside)>1
        heplot1 = (hepred(:,counter+p*16,1));
        heplot2 = (hepred(:,counter+p*16,2));
        plot(concentrations*100, flip(heplot1),'.', 'Color',colorblind(8,:), 'MarkerSize',12)      
        hold on
        if postprocess == 0
            plot((concentrations)*100, (heplot2) ,'.', 'Color',colorblind(1,:), 'MarkerSize',12)
            hold on
        end 
    else
        heplot1 = (hepred(:,counter+p*16));
        plot(concentrations*100, heplot1,'.', 'Color',colorblind(8,:), 'MarkerSize',12)
        hold on 
    end 
    
    %extract experimental predictions 
    
    
    %extract unifac predictions 
    load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
    [~,mixunind] = ismember(mixpred,mixture', 'rows');
    if mixunind ~=0
        heuniplot = he(:,mixunind);
        plot(conc_interval*100, (heuniplot),'k--')
    else 
        [~, mixunind] = ismember([mixpred(3:4) mixpred(1:2)],mixture','rows');
        if mixunind ~=0
            heuniplot = he(:,mixunind);
            plot(conc_interval*100, (heuniplot),'k--')
        else
            disp(mixpred)
            disp('Not found')
        end 
    end 
    %create the title 
    [~,index1]= ismember(mixpred(:,[1,2]),codesNames,'rows');
    [~,index2]= ismember(mixpred(:,[3,4]),codesNames,'rows');
    compoundname1 = compoundnames{index1};
    compoundname2 = compoundnames{index2};
    title(strcat(compoundname1, " & ", compoundname2))% work on making this words - dictionary of compounds and mixture numeric? 
    %title(strcat(label_func(find(strcmp(func_group, {(mixpred(1,1))}))), ' ', num2str(mixpred(1,2)), '+', label_func(find(strcmp(func_group, {num2str(mixpred(1,3))}))), ' ',num2str(mixpred(1,4))))) 
    
    hold off
end 
xlabel('Composition of compound 1 (%)')
ylabel('Excess enthalpy (J/mol)')
if length(whichside)>1 && postprocess ==0
    legend('Experimental','Predictions left', 'Predictions right','UNIFAC (Do)', 'Location', 'southeast')
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
%
load(filename, 'comps', 'mixtureT')

func_group = {1,2,3:5, [6,7,11], 8:9, 10, 12:17,18:20, 21, 22};
label_func = {'Alkane', 'Primary alcohol', 'Other alcohol','Cycloalkanes', 'Ketone', 'Alkene', 'Ester', 'Amine','Acid','Aldehyde'};

err_funcgroup = cell(length(func_group),length(func_group));
ard_funcgroup = cell(length(func_group),length(func_group));
smse_funcgroup = zeros(length(func_group),length(func_group));

wsmse_funcgroup = zeros(length(func_group),length(func_group));
AARD_funcgroup = zeros(length(func_group),length(func_group));
no_mix = zeros(length(func_group),length(func_group));

%extracting the correct mixtures, if necessary 

mixtureT = mixtures(filled_ind,:);

for func = 1:length(func_group)
    for func2 = 1:length(func_group)
        %extract arrays of groups 
        funcgroup1 = func_group{func};
        funcgroup2 = func_group{func2};
        %find mixtures with this combination
        indexloop = 1;
        index=[];
        for i =1:length(funcgroup1)
            indextemp = find(mixtureT(:,1)==funcgroup1(i));
            if ~isempty(indextemp) %&& index~=0
                index(indexloop:length(indextemp)+indexloop-1,1) = indextemp;
                indexloop = length(index)+indexloop;
            end
%             clear indextemp 
%             indextemp = find(mixture(:,3) == funcgroup1(i));
%             if ~isempty(indextemp) %&& index2~=0
%                 index(indexloop:length(indextemp)+indexloop-1,1) = indextemp;
%                 indexloop = length(index)+indexloop;
%             end
        end 
        %index = (index1,index2);
        indexloop = 1;
        index2=[];
        for i =1:length(funcgroup2)
            indextemp = find(mixtureT(:,3)==funcgroup2(i));
            if ~isempty(indextemp) %&& index2~=0
                index2(indexloop:length(indextemp)+indexloop-1,1) = indextemp;
                indexloop = length(index2)+indexloop;
            end
%             clear indextemp 
%             indextemp = find(mixture(:,1) == funcgroup2(i));
%             if ~isempty(indextemp) %&& index2~=0
%                 index2(indexloop:length(indextemp)+indexloop-1,1) = indextemp;
%                 indexloop = length(index2)+indexloop;
%             end
        end 

        indices = intersect(index,index2);
        indices(indices==0) =[];
        %export errors and AARDs
        if length(indices)>1
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
end 

% Heatmap of the errors per functional group 

SMSEplot = smse_funcgroup;
AARDplot=AARD_funcgroup.*100;
tempvar = triu(SMSEplot);
tempvar(tempvar==0)=nan;
%create a table 
[row,col] = find(~isnan(tempvar));
clear smse_tbl aard_tbl wsmsetbl funcgroup1tbl funcgroup2tbl
for counter = 1:length(row)
    funcgroup1tbl{counter} = label_func{row(counter)};
    funcgroup2tbl{counter} = label_func{col(counter)};
    smse_tbl(counter) = SMSEplot(row(counter), col(counter))+SMSEplot(col(counter),row(counter));
    aard_tbl(counter) = AARDplot(row(counter),col(counter))+AARDplot(col(counter),row(counter));
    wsmsetbl(counter) = wsmse_funcgroup(row(counter),col(counter)) + wsmse_funcgroup(col(counter),row(counter));
end 

% heatmap 
tbl = table(funcgroup1tbl', funcgroup2tbl',smse_tbl',aard_tbl',wsmsetbl', 'VariableNames', ["Functional group 1", "Functional group 2", "SMSE (J/mol)","AARD (%)","wSMSE (J/mol)" ]);
% heatmap 
f1 = figure(1);
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'SMSE (J/mol)', 'ColorMethod', 'none');
% h.YDisplayData = flip(label_func);
% h.XDisplayData = label_func;
h.Colormap=winter;
h.FontSize = 14;

% heatmap 
f2 = figure(2);
clf
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'AARD (%)', 'ColorMethod', 'none');
% h.YDisplayData = flip(label_func);
% h.XDisplayData = label_func;
h.Colormap=winter;
h.FontSize = 14;
%% Parity plot and error distribution plots 
%Import all the UNIFAC predictions 
plotuni = 1;
load('colorblind_colormap.mat')
load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
for i = 1:length(filled_ind)
    mixpred = mixtureT(i,:);
    [~,mixunind] = ismember(mixpred,mixture', 'rows');
    if mixunind ~=0
        heuni(:,i) = he(5:5:96,mixunind); % check these indices 
        indexuni(i)=1;
    else
        [~,mixunind] = ismember([mixpred(3:4) mixpred(1:2)],mixture', 'rows');
        if mixunind ~=0
            heuni(:,i) = he(5:5:96,mixunind); % check these indices 
            indexuni(i)=1;
        else 
            indexuni(i) = 0;
        end 
        
    end 
    
end 
 

%parity plot 
f4 = figure(4);
f4.Position=[10 10 500 400];
clf
%choose between plot and loglog 
plotx =Truth(:,indicespred,whichside);
ploty = hepred(:,indicespred,whichside);
plot(plotx(:),ploty(:),'.','Color', 'b', 'MarkerSize',6)
hold on 
if plotuni
    if length(whichside) >1
        plot(plotx(:), [heuni(:); heuni(:)],'.','Color', 'k', 'MarkerSize',6)
    else 
        plot(plotx(:), [heuni(:)],'.','Color', 'k', 'MarkerSize',6)
    end 
        
end 
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

% error distribution 2-way
f5 = figure(5);
f5.Position = [10 20 400 300];
clf
%choose the number of bins 
%ploterrhist(errorHE(:), 'bins', 20)
histogram(errorHE(:), 20, 'FaceColor', colorblind(6,:), 'BinLimits', [-1000 1000])
ylabel('Instances','FontSize',13, 'Color', 'k')
xlabel('MCM error (J/mol)', 'FontSize',13, 'Color', 'k')
% 3way histogram 
f6 = figure(6);
f6.Position = [10 20 450 300];
clf
edges = {-1000:100:1000; 0:5:100}; % bin edges
data = [errorHE(:) abs(aardHE(:)).*100]; 
hist3(data, 'Edges', edges,'CdataMode','auto','FaceColor','interp')
zlabel('Instances','FontSize',13, 'Color', 'k')
xlabel('Error (J/mol)', 'FontSize',13, 'Color', 'k')
ylabel('ARD (%)' ,'FontSize', 13, 'Color', 'k')
colormap(colorblind(flip([9,1,8,6,3]),:))

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
xlabel('Compound 1')
ylabel('Compound 2')
zlabel('Composition of compound 1 (%)', 'FontSize',12)
indices = 1:4:length(compoundnames);
xticks(indices)
yticks(indices)
xticklabels(compoundnames(indices))
yticklabels(compoundnames(indices))
%% Classification of correct predictions per system 
% colours of data points will show which one is better 
%started 04/11/2022 - check if it works 
%scatter(x,y,pointsize,color)
fnind = 1; %can change this here 
load('colorblind_colormap.mat')
errorHEsystem1 = reshape(errorHEsystem{fnind}, length(filled_ind),2);
%errorUNIsystem = reshape(errorUNIsystem, length(filled_ind),2);
leftcorrect = (errorHEsystem1(:,1)<errorHEsystem1(:,2));% logical output 
errorHEbest = errorHEsystem1(:,1).*(leftcorrect) + errorHEsystem1(:,2).*(1-leftcorrect);

unifaccorrect = (errorHEbest>errorUNIsystem'); %logical output 
numleftcorrect = sum(find(leftcorrect));
numunicorrect = sum(find(unifaccorrect));
%figure out x and y later 
leftC = nan(size(X,[1,2]));
uniC = nan(size(X,[1,2]));
for i = 1:length(filled_ind)
    leftC(row(i),col(i)) = leftcorrect(i);
    uniC(row(i),col(i)) = unifaccorrect(i);
end 
figure(6)
clf
[x,y] = meshgrid(1:size(X,1), 1:size(X,2));
pcolor(x,y,leftC)
xlabel('Compound 1', 'FontSize',12)
ylabel('Compound 2', 'FontSize',12)
c = colorbar('Ticks', [0,1], 'TickLabels', ['Upper'; 'Lower']); %upper triangular = left, 0 is false
colormap([colorblind(6,:); colorblind(8,:)])
indices = 1:2:length(compoundnames);
xticks(indices)
yticks(indices)
xticklabels(compoundnames(indices))
yticklabels(compoundnames(indices))
figure(7)
clf
[x,y] = meshgrid(1:size(X,1), 1:size(X,2));
pcolor(x,y,uniC)
xlabel('Compound 1', 'FontSize',12)
ylabel('Compound 2', 'FontSize',12)
c = colorbar('Ticks', [0,1], 'TickLabels', ["UNIFAC"; "MCM"]); % 0 is false
colormap([colorblind(6,:); colorblind(8,:)])
indices = 1:2:length(compoundnames);
xticks(indices)
yticks(indices)
xticklabels(compoundnames(indices))
yticklabels(compoundnames(indices))
%% Compounds and numeric matches
[compoundnames, codesNames] = findnames([mixturepred(:,[1,2]); mixturepred(:,[3,4])]);

%% functions 
function [compoundnames,codesout] = findnames(codesin)
    compounds = {'Methane','Ethane', 'Propane', 'Butane', 'Pentane', 'Hexane', 'Heptane', 'Octane', 'Nonane', 'Decane', 'Dodecane', 'Methanol', 'Ethanol', '1-Propanol', '1-Butanol', '1-Pentanol', '1-Hexanol', '1-Heptanol', '1-Octanol', '1-Nonanol', '1-Decanol', '2-propanol', '2-butanol', '2-pentanol', '2-hexanol','2-octanol', 'Isobutanol', 'Tertbutylalcohol', 'Cyclopentane', 'Cyclohexane','Cycloheptane', 'Cyclooctane','Benzene', 'Toluene', 'Ethanal', 'Propanal', 'Butanal', 'Pentanal', 'Hexanal', 'Ethene', '1-Pentene', '1-Hexene', '1-Heptene', '1-Octene', '1-Nonene', '1-Decene', 'Acetone', 'Butanone', '2-Propanone', '2-Hexanone', '2-Heptanone', '2-Octanone','Formic acid', 'Acetic acid', 'Propionic acid', 'Butyric acid', 'Pentanoic acid', 'Hexanoic acid', 'Heptanoic acid', 'Octanoic acid', 'Nonanoic acid', 'Decanoic acid', '3-Pentanone', '3-Hexanone', '3-Heptanone', '3-Octanone', 'Methyl formate', 'Methyl acetate', 'Methyl propionate', 'Methyl butyrate', 'Methyl pentanoate', 'Methyl hexanoate', 'Methyl benzoate', 'Ethyl benzoate', 'Ethyl formate', 'Ethyl acetate', 'Ethyl propionate', 'Ethyl butyrate', 'Ethyl pentanoate', 'Ethyl hexanoate', 'Ethyl heptanoate', 'Ethyl octanoate' , 'Propyl formate','Propyl acetate', 'Propyl propionate', 'Propyl butyrate', 'Butyl formate', 'Butyl acetate', 'Butyl butyrate', 'Pentyl acetate', 'Propylamine', 'Butylamine', 'Pentyalamine', 'Hexylamine', 'Heptylamine', 'Octylamine', 'Nonylamine', 'Decylamine', 'Aniline', 'Benzylamine'};
    codes = [ones(11,1), ([1:10,12])'; 2*ones(10,1), (1:10)'; 3*ones(5,1), ([3:6,8])'; 4, 4; 5, 4;11*ones(4,1), (5:8)'; 6, 0; 7, 0; 22*ones(5,1), (2:6)'; 10*ones(7,1), ([2,5:10])'; 8*ones(6,1), (3:8)'; 21*ones(10,1), (1:10)'; 9*ones(4,1), (5:8)';12*ones(6,1), (1:6)'; 17*ones(2,1), (1:2)'; 13*ones(8,1), (1:8)'; 14*ones(4,1), (1:4)'; 15*ones(3,1), ([1,2,4])'; 16, 1; 18*ones(8,1), (3:10)'; 19, 0; 20, 0];
    indices = find(ismember(codes, codesin,'rows'));
    for i = 1:length(indices)
        compoundnames{i} = compounds{indices(i)};
        codesout(i,:) =codes(indices(i),:);
    end 
end 

function [wmse]=find_wmse_error(errors, count)
    errors = reshape(errors,prod(size(errors)),1);
    perc5=prctile(errors,5,'all');
    perc95=prctile(errors,95,'all');
    errors(errors<perc5)=perc5;
    errors(errors>perc95)=perc95;
    wmse = (sum((errors).^2))/count;
end 

function [heprednew,uncertainty] = postprocesshe(hepred,concentrations)
    heprednew = zeros(size(hepred,1),size(hepred,2));
    
    uncertainty = zeros(length(size(hepred,2)),length(concentrations));
    threshold = 0.8*(size(hepred,1)*2-1)/sqrt(size(hepred,1)*2);

    for index = 1:(size(hepred,2))
        temppred = hepred(:,index,:);
        if size(hepred,3)>1
            xvar = [concentrations,concentrations];
        else 
            xvar = concentrations;
        end 
        temppred = temppred(:);
        %remove abnormal values and zeros 
        temppred(abs(temppred)>1e10)=nan;
        temppred(abs(temppred)<1e-3)=nan;
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
        [p,S,mu] = polyfit(xvar, temppred,6);
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