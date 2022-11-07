%% Francesca Middleton 27/09/2022

%discovering the percent correct prediction of the sign predictions 

clc
clear
Temperature=298.15;
concen = 1;
ways = 2;
fillmethod = 'dia';
parafac = 'INDAFAC';
%if only one side of the 2-way array is considered or both (1:2)
whichside =1:2;
postprocess = 0;

%% Xsign
%load('2wayINDAFAC-LOOCV-Xsign-maxiter=20000-T=298.15-c=19-fillmethod=dia-21-Oct-2022.mat')

%load('2wayINDAFACthresholded-LOOCV-Xsign-maxiter=20000-T=298.15-c=19-fillmethod=dia-21-Oct-2022.mat')
%load the sign results
%load('2wayINDAFAC-25comps-LOOCV-Xsign-maxiter=20000-T=298.15-c=19-fillmethod=dia-18-Oct-2022.mat')
%load('2wayINDAFAC-All-LOOCV-Xsign-maxiter=20000-T=298.15-c=19-fillmethod=dia-17-Oct-2022.mat')
%change the date manually 
%filenameimport = strcat(num2str(ways),'way',(parafac),'-All-LOOCV-Xsign-maxiter=20000-T=',num2str(Temperature),'-c=',num2str(concen),'-fillmethod=',fillmethod,'-',date,'-2022.mat');
%load(filenameimport)
%Tablesign = {filenameimport; filename; fillmethod; num2str(orth)};
%writetable(cell2table(Tablesign),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'U2', 'WriteVariableNames',false, 'AutoFitWidth', false)
% signPrediction = sign(Xm_boot);
% signPrediction(signPrediction ==0)=nan;
% indicesSign = find(~isnan(signPrediction(1,1,:)));
% signPrediction = signPrediction(:,:,indicesSign);
% Xtemp = Xsign;
% Xtemp(X==0)=nan;
% filled_ind = find(~isnan(Xtemp(:,:,1))); % not equal to zero and not nan
% filled_indices = find(~isnan(X(:,:,1)));
% [row,col] = find(~isnan(Xtemp(:,:,1)));
% 
% errors = zeros(length(fns),length(concentrations),length(filled_ind)/2,2);
% signPreds = zeros(length(fns),length(concentrations),length(filled_ind)/2,2);
% correct = zeros(length(fns),length(concentrations));
% %find errors and how many were correctly predicted 
% 
% Preds = [Xm_boot(1,1,1:length(filled_ind))];%; Xm_boot2(1,1,:)];
% %Preds = reshape(Preds,2,[]);
% indiceskeep = find(Preds(1,:)); % not equal to zero 
% for fn = 1:length(fns)
%     for c = 1:length(concentrations)
% 
%         if length(whichside)>1
%             Preds(1,:) = Xm_boot(c,fn,1:length(filled_ind));
%             Preds(2,:) = Xm_boot2(c,fn,1:length(filled_ind));
%         elseif whichside ==1
%             Preds(1,:) = Xm_boot(c,fn,1:length(filled_ind));
%         else 
%             Preds(1,:) = Xm_boot2(c,fn,1:length(filled_ind));
%         end 
%         
%         Preds = sign(Preds(:, indiceskeep));
%         signPreds(fn,c,:,whichside) =reshape(Preds,size(signPreds(1,1,:,whichside))); 
%         
%         if length(whichside)>1
%             %Truth: lower half (preds1)
%             tempX = tril(Xscale(:,:,c),-1)+triu(nan(size(Xscale(:,:,c))));
%             filled_ind = find(~isnan(tempX));
%             Truthtemp = Xsign(:,:,c);
% 
%             Truthsign(c,1,:) = Truthtemp(filled_ind);
%             %top half truth
%             tempX = (triu(Xsign(:,:,c),1)+tril(nan(size(Xsign(:,:,c)))))';
%             filled_ind = find(~isnan(tempX));
%             Truthsign(c,2,:) = (tempX(filled_ind));
%         elseif whichside==1
%             %Truth: lower half (preds1)
%             tempX = tril(Xscale(:,:,c),-1)+triu(nan(size(Xscale(:,:,c))));
%             filled_ind = find(~isnan(tempX));
%             Truthtemp = Xsign(:,:,c);
% 
%             Truthsign(c,1,:) = Truthtemp(filled_ind);
%         elseif whichside==2
%             %top half truth
%             tempX = (triu(Xsign(:,:,c),1)+tril(nan(size(Xsign(:,:,c)))))';
%             filled_ind = find(~isnan(tempX));
%             Truthsign(c,2,:) = (tempX(filled_ind));
%         end 
%         
%         errors(fn,c,:,whichside) =reshape(reshape(Truthsign(c,whichside,:),size(Preds))-Preds, size(errors(1,1,:,whichside)));
%         correct(fn,c)= sum(sum(errors(fn,c,:,whichside)==0))/length(Preds)/(length(whichside))*100;
%         
%     end 
% end 
%export the results 
%writetable(array2table(fns'),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'T11', 'WriteVariableNames',false, 'AutoFitWidth', false)
%writetable(array2table(correct),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'U11', 'WriteVariableNames',false, 'AutoFitWidth', false)
%% If no sign predictions exist 
load('3waPARAFAC-Small-concslices-LOOCV-Xscale-maxiter=20000-Temp=298.15-fillmethod=reg-04-Nov-2022.mat')
%load('2waySVD-50comps-scale1-threshold60par-moreiter-huberloss1-LOOCV-Xscale-T=298.15-fillmethod=rec-02-Nov-2022.mat')
%load('2waySVD-50comps-scale0-threshold60par-moreiter-huberloss1-LOOCV-Xscale-T=298.15-fillmethod=reb-01-Nov-2022.mat')
%load('2waySVD-Small-threshold60par-moreiter-huberloss1-LOOCV-Xscale-T=298.15-fillmethod=reg-01-Nov-2022.mat')
%load('2waySVD-Small-threshold60par-moreiter-huberloss1-LOOCV-Xscale-T=298.15-fillmethod=avr-01-Nov-2022.mat')
%load('2waySVD-Small-threshold40par-huberloss1-LOOCV-Xscale-T=298.15-fillmethod=avr-01-Nov-2022.mat')
%load('2waySVD-Small-thresholdpar-huberloss1-LOOCV-Xscale-T=298.15-fillmethod=avg-31-Oct-2022.mat')
%load('2wayINDAFAC-Small-threshold50parallel-LOOCV-Xscale-T=298.15-fillmethod=dia-28-Oct-2022.mat')
%load('2wayINDAFAC-Small-threshold50parallel-LOOCV-Xscale-T=298.15-fillmethod=uni-27-Oct-2022.mat')
%load('2waySVD-26comps-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=avg-26-Oct-2022.mat')
%load('2waySVD-Half-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-25-Oct-2022.mat')
%load('2waySVD-26comps-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=dia-25-Oct-2022.mat')
%load('2waySVD-Half-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-25-Oct-2022.mat')
%load('2wayPARAFAC-LOOCV-Half-Xscale-maxiter=20000-T=298.15-c=6-fillmethod=uni-25-Oct-2022.mat')
 %load('2waySVD-Alll-LOOCV-Xscale-maxiter=20000-T=298.15-c=1-fillmethod=dia-23-Oct-2022.mat')
%insert scale results filename here 
%load('2wayINDAFACcentred-All-LOOCV-Xscale-maxiter=20000-T=298.15-c=19-fillmethod=dia-17-Oct-2022.mat')
%load('2wayINDAFACcentred-All-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-18-Oct-2022.mat')
%load('2wayPARAFACthresholded-LOOCV-Xscale-maxiter=20000-T=298.15-c=19-fillmethod=dia-24-Oct-2022.mat')
%load('2waySVDscaled+centred+thresholdZ-All-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-23-Oct-2022.mat')
%load('2waySVDscaled+centred-All-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-20-Oct-2022.mat')
%load('2wayINDAFAC-All-LOOCV-Xscale-maxiter=2500-T=298.15-c=19-fillmethod=uni-20-Oct-2022.mat')
%load('2wayINDAFAC-All-LOOCV-Xorig-maxiter=20000-T=298.15-c=9-fillmethod=uni-13-Oct-2022.mat')
%finds the true sign values
X(X==0)=nan;
% for counter =1:size(X,3)
%     Xnew(:,:,counter,Tempind) = remove_nan2(X(:,:,counter,Tempind));
% end 

%X = Xnew;
Xsign = sign(X);
concentrations = conc_interval;
tempX = tril(X(:,:,1,Tempind),-1)+triu(nan(size(X(:,:,1,Tempind))));
[row,col] = find(~isnan(tempX));
filled_ind =  find(~isnan(tempX));
         

for index = 1:length(filled_ind)
    Truthsign(:,index,1) = Xsign(row(index),col(index),:,Tempind);
    Truthsign(:,index,2) = (Xsign(col(index),row(index),:,Tempind)); 
end 



%% Xscale
date ='13-Oct';

%load('2waySVD-Half-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-25-Oct-2022.mat')
%load('2waySVD-26comps-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=dia-13-Oct.mat')
%clear Truth Truthtemp scalePreds hepred scalePrediction mse wmse smse wsmse smse
%load('2waySVD-Half-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-25-Oct-2022.mat')
%load('2wayINDAFACthreshold-All-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-21-Oct-2022.mat')
%load('2wayPARAFAC-LOOCV-Half-Xscale-maxiter=20000-T=298.15-c=6-fillmethod=uni-25-Oct-2022.mat')
%load('2waySVD-Half-LOOCV-Xscale-maxiter=25000-T=298.15-c=18-fillmethod=dia-13-Oct-2022.mat')
%load('2wayPARAFACthresholded-LOOCV-Xscale-maxiter=20000-T=298.15-c=19-fillmethod=dia-24-Oct-2022.mat')
%load('2waySVDscaled+centred+thresholdZ-All-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-23-Oct-2022.mat')
%load('2waySVDscaled+centred-All-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-20-Oct-2022.mat')
%load('2wayINDAFAC-25CompsMoreFns-LOOCV-Xscale-maxiter=20000-T=298.15-c=19-fillmethod=dia-19-Oct-2022.mat')
%load('2waySVD-Alll-LOOCV-Xscale-maxiter=20000-T=298.15-c=9-fillmethod=dia-17-Oct-2022.mat')
%load('2waySVD-Alll-LOOCV-Xscale-maxiter=20000-T=298.15-c=14-fillmethod=dia-17-Oct-2022.mat')
%load('2wayINDAFAC-All-LOOCV-Xscale-maxiter=2500-T=298.15-c=19-fillmethod=uni-20-Oct-2022.mat')
%load('2wayINDAFACthresholded-LOOCV-Xscale-maxiter=20000-T=298.15-c=19-fillmethod=dia-21-Oct-2022.mat')
%load('2wayINDAFACcentred-All-LOOCV-Xscale-maxiter=20000-T=298.15-c=19-fillmethod=dia-17-Oct-2022.mat')
%load('2wayINDAFACcentred-All-LOOCV-Xscale-maxiter=25000-T=298.15-c=19-fillmethod=uni-18-Oct-2022.mat')
%load('2wayINDAFAC-25comps-LOOCV-Xscale-maxiter=20000-T=298.15-c=19-fillmethod=dia-18-Oct-2022.mat')
% filenameimport = strcat(num2str(ways),'way',(parafac),'centre-All-LOOCV-Xscale-maxiter=20000-T=',num2str(Temperature),'-c=',num2str(concen),'-fillmethod=',fillmethod,'-',date,'-2022.mat');
% load(filenameimport)
%Tablescale = {filenameimport; filename; fillmethod; num2str(orth)};
%writetable(cell2table(Tablescale),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'B2', 'WriteVariableNames',false, 'AutoFitWidth', false)
scalePrediction = Xm_boot;
scalePrediction(scalePrediction ==0)=nan;
indicesScale = find(~isnan(scalePrediction(1,1,:)));
%indicesInt = intersect(indicesSign,indicesScale); % k that you can use to judge He predictions 
Xtemp = Xs(:,:,1,Tempind);
tempX = tril(Xtemp,-1)+triu(nan(size(Xtemp)));


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
                
                Truthscale(:,index,2) = Xscale(col(index),row(index),:);
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

%export 
% writetable(array2table(fns'),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'A10', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(fns'),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'A25', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(fns'),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'A40', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(mse),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'B10', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(wmse),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'B25', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(aard),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'B40', 'WriteVariableNames',false, 'AutoFitWidth', false)
%% He preds
fnind=0;
load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
[~,mixunind] = find(ismember(mixture',mixtures(filled_ind,:), 'rows')); % hopefully this extracts what we want 
heunifacpred = he(5:5:96, mixunind);
mixtureuni = mixture;
for fnind =6%1:length(fns)

    rSign = 1;
    rScale = fnind;
    
    
    %hepred = signPreds(rSign,:,:,:).*exp(signPreds(rSign,:,:,:).*scalePreds(rScale,:,:,:));
    
    hepred = Truthsign(:,:,whichside).*exp(reshape(scalePreds(rScale,:,:,whichside),size(Truthsign(:,:,whichside))).*Truthsign(:,:,whichside));
    errorUNI(:,:,whichside) = hepred - heunifacpred; 
    errorUNIsystem = mean(errorUNI,2);
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
    indicespred = find(~isnan(hepred(concen,:,1)));
    if whichside ==2
        errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred));
    else 
        errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred,whichside));
    end 
    errorHE(isnan(errorHE)) = 0;
    errorHEsystem = mean(errorHE,2);
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
%find the minimum for each concentration
% for c =1:length(concentrations)
%     ind = find(min(wsmsec(:,c))==wsmsec(:,c));
%     if isempty(ind)
%         ind = fns(fnind);
%     end 
%     tableout(:,c) = [rSign; fns(ind); smsec(ind,c); aardc(ind,c); wsmsec(ind,c); reshape(percless(ind,c,:),5,1)];
% end 
%overall errors 
% writetable(array2table(fns(fn)'),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'AL6', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(wsmse),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'AM6', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(wsmsec),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'AN6', 'WriteVariableNames',false, 'AutoFitWidth', false)
% writetable(array2table(tableout),'2waycompletionScaled.xlsx', 'Sheet', num2str(Temperature),'Range', 'AM24', 'WriteVariableNames',false, 'AutoFitWidth', false)
%%
figure(1)
clf
index1plot=1;
semilogy(fns(index1plot:end),smse(index1plot:end),'.', 'MarkerSize',14)
hold on 
semilogy(fns(index1plot:end),wsmse(index1plot:end),'.', 'MarkerSize',14)

ylabel('Error (J/mol)')
xlabel('Number of factors')
legend('SMSE','wSMSE')
%%
clf
index1plot=1;
%choose mse or wmse
semilogy(fns(index1plot:end),smse(index1plot:end),'.','MarkerSize',12)

%semilogy(fns(index1plot:end),wsmse(index1plot:end),'.','MarkerSize',12)

ylabel('SMSE (J/mol)')

xlabel('Number of factors')

%% Define mixture you wish to predict
figure(2)
clf

p = 0;


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
    plot(conc_exp*100, heexp,'r.', 'MarkerSize',12)
    hold on
    
    %predictions made by the model 
    if length(whichside)>1
        heplot1 = (hepred(:,counter+p*16,1));
        heplot2 = (hepred(:,counter+p*16,2));
        plot(concentrations*100, flip(heplot1),'b.', 'MarkerSize',12)      
        hold on
        if postprocess == 0
            plot(concentrations*100, (heplot2) ,'c.', 'MarkerSize',12)
            hold on
        end 
    else
        heplot1 = (hepred(:,counter+p*16));
        plot(concentrations*100, heplot1,'b.', 'MarkerSize',12)
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
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'wSMSE (J/mol)', 'ColorMethod', 'none');
h.Colormap=winter;
h.FontSize = 14;

%% Parity plot and error distribution plots 
%Import all the UNIFAC predictions 
load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
for i = 1:length(filled_ind)
    mixpred = mixtureT(i,:);
    [~,mixunind] = ismember(mixpred,mixture', 'rows');
    if mixunind ~=0
        heuni(:,i) = he(5:5:96,mixunind); % check these indices 
    end 
end 
figure(4) 
clf 

%parity plot 
subplot(2,1,1)
%choose between plot and loglog 
loglog(Truth(:,indicespred,whichside),hepred(:,indicespred,whichside), 'bo', 'MarkerSize',12)
hold on 
loglog(Truth(:,indicespred,1), heuni, 'ko', 'MarkerSize',12)
loglog([0,max(max(max(Truth)))], [0, max(max(max(Truth)))], 'Color',[.7 .7 .7])

hold off 
xlabel('Predictions (J/mol)')
ylabel('Experimental (J/mol)')

%error distribution 
subplot(2,1,2)
%choose the number of bins 
ploterrhist(errorHE, 'bins', 20)

%% Classification of correct predictions per system 
% colours of data points will show which one is better 
%started 04/11/2022 - check if it works 
%scatter(x,y,pointsize,color)

leftcorrect = errorHEsystem(:,1)<errorHEsystem(:,2);% logical output 
errorHEbest = errorHEsystem(leftcorrect,1) + errorHEsystem(1-leftcorrect,2);
unifaccorrect = errorHEbest>errorUNIsystem; %logical output 
numleftcorrect = sum(find(leftcorrect));
numunicorrect = sum(find(unifaccorrect));
%figure out x and y later 
scatter(1:length(filled_ind),1:length(filled_ind),[],leftcorrect, 'filled')
scatter(1:length(filled_ind),1:length(filled_ind),[],unifaccorrect, 'filled')

%Classification by machine learning method 
%% Classification of correct predictions per functional group 
%fix all of this before using - 04/11/2022

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
        unifac_funcgroup{func,func_2} = sum(unifaccorrect(indices));
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
    threshold = 0.6*(size(hepred,1)*2-1)/sqrt(size(hepred,1)*2);

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
