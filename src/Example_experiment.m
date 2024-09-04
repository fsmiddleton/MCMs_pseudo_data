%% 2-way arrays completed in parallel using the functions in this library
% Francesca Middleton, 2024-09-4


% Parameters
r = [3:10]; % ranks
T = 298.15; % Temperature (K)
fillmethod = 'uni'; % Filling metod used
maxiter = 50000; % Maximum iterations used
thresholdperc = 0.5; % Vlaue of the coherence constraint 
filename = strcat('HEData3wayPolyAll',num2str(T), '.mat'); % Filenmae constructuted for input 
load(filename)

%% Perform parallel matrix completion on the array
[filenamesave,Xm_boot,Xm_boot2,r, X, X,conc_interval,filename, filled_ind] = completion_2way_par(r,fillmethod,maxiter,filename,thresholdperc);
save(filenamesave)

%% Perform some analysis on the results  
% Load the original data 
load(filename)
% Find filled indices in the original array 
Xtemp(X==0)=nan;
tempX = Xtemp(:,:,1);
indend = 00; %this is valid for when compounds were added
tempX(1:indend,1:indend) = nan;
tempX = tril(tempX,-1)+triu(nan(size(tempX)));
[row,col] = find(~isnan(tempX));
filled_ind = find(~isnan(tempX));

% Load the results 
load(filenamesave)
% initiate the array holding predictions 
Preds = [X(1,1:length(filled_ind))];%; Xm_boot2(1,1,:)];
errors = zeros(length(r),size(X,3),length(filled_ind),2);
wsmse = zeros(length(r),size(X,3));
fnind = 0;
% Consider all ranks 
for fn = 1:length(r)
    fnind = fnind+1;
    for c = 1:size(X,3)
        % Predictions 
        Preds(1,:) = Xm_boot(c,fn,1:length(filled_ind));
        Preds(2,:) = Xm_boot2(c,fn,1:length(filled_ind));
    
        % Truth: lower half 
        Truthtemp = X(:,:,c);
        Truth(fnind,c,:,1) = Truthtemp(filled_ind);
        % top half 
        tempX = (triu(X(:,:,c),1)+tril(nan(size(X(:,:,c)))))';
        filled_ind = find(~isnan(tempX));
        Truth(fnind,c,:,2) = (tempX(filled_ind));
        
        errors(fnind,c,:,:) =reshape(reshape(Truth(fnind,c,:,:),size(Preds))-Preds, size(errors(1,1,:,:)));
        wsmse(fnind,c) = find_wmse_error(errors(fnind,c,:,:));
    end 
end 

%% Plot these 

% End