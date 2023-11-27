%% 2-way arrays 
% Francesca Middleton, 2022-03-02

%% Import or create data matrix for testing 
clc
clear
% Import excess enthalpy data for one matrix 
filename = 'TestHEMatrixSMALL15June298.15.xlsx';
%filename = 'TestHEMatrixSMALLreduced15June2022298.15.xlsx';
T1 = readtable(filename, 'Sheet', '0.1'); 
X = table2array(T1);%only use T1, small enough to manage
dim = size(X);
dim1=dim(1);
dim2=dim(2);

% Analysis of the 2-way array 
percmiss = length(find(isnan(X)))/(dim1*dim2)*100;
percobs = length(find(~isnan(X)))/(dim1*dim2)*100;
%parpool
b1 = table2array(readtable(filename, 'Sheet', 'B1'));
b2 = table2array(readtable(filename, 'Sheet', 'B2'));
[~,Locb] = ismember(b2,b1,'rows');
include2 = find(Locb==0);
% all possible components in this matrix 
comps = [b1; b2(include2,:)];

mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end 

%% Plot the matrix generated
clf
%[Xs,missing_ind,filled_linear_ind]=fill_matrix(X,0.5);%missing data 
hm = HeatMap(X);
addXLabel(hm,'Component 1','FontSize',12);
addYLabel(hm,'Component 2','FontSize',12);
view(hm)
histogram(X)
xlabel('Value in the matrix X')
ylabel('Frequency')
%% Rank determined from SVs

% Import other data for most missing to allow variables to be declared
concentrations = 0.1:0.1:0.9;
fn = min([dim1,dim2]);

y1 = zeros(fn,length(concentrations));
y2 = zeros(fn,length(concentrations));
Sts = zeros(fn,fn,length(concentrations));
%necessary for plots
indexplot=20;

i=0;
for c = concentrations
    i=i+1;
    %Import data 
    T = readtable(filename, 'Sheet', num2str(c));
    Xs=table2array(T);
    Xf = filldata2(Xs,'uni',mixtures,c);
    %SVD
    % [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd(X,fn,center,scale,conv,max_iter, use_missing)
    [~,D,~,St,~, ~]=missing_svd(Xf,fn,1,0,1e-8,1,2);
    %Plot
    y1(:,i) = diag(D);
    y2(:,i) = cumsum(y1(:,i))/sum(y1(:,i));
    Sts(:,:,i)=St;
    subplot(2,1,1)
    plot(1:indexplot, y1(1:indexplot,i))
    subplot(2,1,2)
    plot(1:indexplot,y2(1:indexplot,i))
    hold on
end 

hold off 

%% Preprocess the matrices

concentrations = 0.1:0.1:0.9;
corrcol = zeros(length(concentrations), dim1, dim2);
corrrow = zeros(length(concentrations), dim1, dim2);
covX = zeros(length(concentrations), dim1, dim2);
detX = zeros(length(concentrations), dim1, dim2); % = product of singular values 
trX  = zeros(length(concentrations));
for c = 1:length(concentrations)
    i=i+1;
    %Import data 
    T = readtable(filename, 'Sheet', num2str(c/10));
    Xs=table2array(T);
    X = fill_data(Xs);
    corrcol(c, :, :) = corr(X);
    corrrow(c, :, :) = corr(X');
    covX(c, :, :) = cov(X);
    detX(c) = det(reshape(covX(c,:,:), [dim1,dim2])); % = product of singular values 
    trX(c) = trace(reshape(covX(c,:,:), [dim1,dim2])); % = sum of singular values
end 
 
%heatmap(abs(tril(corrcol)))
%% LOOCV to find rank using mse  
%time the code 
tic
% declare missing % used and the file from which to extract information 
fns =6:9;
concentrations=0.1:0.1:0.9;
Xs=X;
missing_ind = find(isnan(Xs));
% largest length of filled_linear_ind is for 0.3
filled_ind = find(~isnan(Xs));

% choose whether to use mse or wmse to find the optimal rank 
winsorized_mse = 1;
% declare vars for analysis 
mse_LOOCV = zeros(length(concentrations),length(fns));
wmse_LOOCV = zeros(length(concentrations),length(fns));
RAD_LOOCV = zeros(length(concentrations),length(fns)); % relative absolute deviation 
min_mse_LOOCV = zeros(length(concentrations),1);
min_wmse_LOOCV = zeros(length(concentrations),1);
RAD_LOOCV_best = zeros(length(concentrations),1);
X_pred_LOOCV = zeros(dim1,dim2,length(filled_ind));
LOOCV_removed_col = zeros(length(filled_ind),1);
LOOCV_removed_row = zeros(length(filled_ind),1);
SVs = zeros(dim1,length(concentrations));

maxiter = 100000;
conv = 1e-10;
Xm_boot=zeros(length(concentrations), length(fns), length(filled_ind));
%concentrations=0.1;

for c = 1:length(concentrations)
    conc = concentrations(c);
    T = readtable(filename, 'Sheet', num2str(conc));
    Xs = table2array(T);
    Xs = remove_nan2(Xs);
    Xfilled = Xs;
    [row,col] = find(~isnan(Xs));
    filled_ind = find(~isnan(Xs));
    disp('Concentration of component 1')
    disp(c)
    %declare vars with size dependent on the array used 
    %loop through ranks
    disp('fn')
    fnind=0;
    for fn = fns 
        disp(fn)
        fnind = fnind + 1; 
        error_LOOCV = zeros(length(filled_ind),1);
        RAD = zeros(length(filled_ind),1);
        parfor k = 1:length(filled_ind') %filled_linear_ind must be a row vector for a for loop    
            filled_index = filled_ind(k);
            % remove a point from Xs
            X_b = Xs;
            X_b(filled_index) = nan;
            if find(~isnan(X_b(:,col(k)))) & find(~isnan(X_b(row(k),:))) & Xs(filled_index)~=0 %ensure at least one value in each column and row
                LOOCV_removed_col(k,c) = col(k);
                LOOCV_removed_row(k,c) = row(k);
                %perform iterative PCA on the slightly more empty matrix 
                %  [S,V,D,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter, usemissing)
                [~,~,~,~,Xpred, ~,iter(k)]=missing_svd(X_b,fn,1,1,conv,maxiter,1);
                X_pred_LOOCV(:,:,k)=Xpred;
                error_LOOCV(k) = Xs(filled_index)-Xpred(filled_index);
                
                Xm_boot(c, fnind, k) = Xpred(filled_index);
                if Xs(filled_index)~=0
                    RAD(k) = error_LOOCV(k)/Xs(filled_index);
                end
                
            end
        end
        % mse for this composition and rank 
        mse_LOOCV(c,fnind)= sum(error_LOOCV.^2)/length(error_LOOCV);
        wmse_LOOCV(c, fnind) = find_wmse_error(error_LOOCV, length(filled_ind'));
        %absolute average deviation
        RAD_LOOCV(c,fnind) = (sum(abs(RAD)))/length(RAD);
    end
    
    % find the optimal rank 
    if winsorized_mse ==1
        fn_LOOCV(c) = find(wmse_LOOCV(c,:) == min(wmse_LOOCV(c,:)));
    else 
        fn_LOOCV(c) = find(mse_LOOCV(c,:) == min(mse_LOOCV(c,:)));
    end
    [U,D,V,St,Xpred, SVsmiss,it]=missing_svd(Xs,fn_LOOCV(c),1,1,conv,maxiter,1);
    y = diag(SVsmiss);
    SVs(:,c)=y;
    min_mse_LOOCV(c) = mse_LOOCV(c,fn_LOOCV(c));
    min_wmse_LOOCV(c) = wmse_LOOCV(c,fn_LOOCV(c));
    RAD_LOOCV_best(c) = RAD_LOOCV(c, fn_LOOCV(c));
    % use the best rank to find error bars 
end 
Results = table(["Concentration of component 1"; (concentrations')],[ "SVs";fn_LOOCV'],[ "MSE";min_mse_LOOCV],[ "wMSE";min_wmse_LOOCV]);
disp(Results)
toc 
% Save variables from the run 
filenametemp = strcat('2wayrunPCAPolySmallLOOCV-',date, '.mat');
save(filenametemp)
%retrieve data using load(filename.mat)
%% Plot errors 
conc=0.1;
c=conc*10;
plot(fns, (RAD_LOOCV(c,:)))

%% Find optimal rank of the model for the filled data 


%% Import .m file with 3-way array
clc
clear
% decide on which interval and temperature to evaluate 
interval = 0.05;
T =308.15;
% choose which array you want to work with 
%filename = 'HEData3wayPolyMid-0.05-298.15.mat';
%filename = 'HEArrayPolySmall298.15.mat';
filename = strcat('HE3wayPolyAll',num2str(T),'.mat');
 load(filename, 'HE_data_sparse',  'comps', 'mixtureT','mixture', 'Temps')
%filename = strcat('HEData3wayPolySmall-',num2str(interval), '-', num2str(T), '.mat');
%filename = strcat('HEData3wayPolyAll-',num2str(interval), '-', num2str(T), '.mat');
%load(filename)
conc_interval = interval:interval:(1-interval);
indend = size(HE_data_sparse,1);
X = HE_data_sparse(1:indend,1:indend,1:length(conc_interval));
Xs = X;
Xsign = sign(X);
Xscale = Xsign.*log(Xsign.*X);
dim = size(X);
dim1 = dim(1);
dim2 = dim(2);
dim3 = dim(3);
percmiss = length(find(isnan(X)))/(dim1*dim2*dim3)*100;
percobs = length(find(~isnan(X)))/(dim1*dim2*dim3)*100;
tempX = tril(X(:,:,1),-1)+triu(nan(size(X(:,:,1))));
[row,col] = find(~isnan(tempX));
filled_ind = find(~isnan(tempX));
% define all the mixtures in the linear indices of the array, based on the
% components in the 3-way array file 
comps = comps(1:indend,:);
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end 
Xfilledall = filldata3(X,'uni', mixtures,conc_interval, 'none', T);

%necessary vars 
dim = size(Xs);
dim1=dim(1);
dim2=dim(2);
missing_ind = find(isnan(Xs));
% largest length of filled_linear_ind is for 0.3
filled_ind = find(~isnan(Xs));
% uses all the available data to find the rank and the first section of code to run 
% remove data or fill a matrix 
%filename specified in first block of code for import ==1 
remove = 0; % 1=remove to create sparse matrix, 2 = fill a sparse matrix, 0 = no change to matrix
remove_ind = 1:(dim1*dim2);
 
%Time to choose how to find the best rank
Gavish =1;
winsorized_mse = 0; %1=use wmse. Used for iterative PCA 
 % Choose sparsities used for finding the rank 
concentrations = 0.05:0.05:0.95;
%ranks to try 
fns = 1:1:10;
%choose whether to reorder matrix or not, 1 = true 
reorder = 0;
        
%intialise the metrics to analyse each sparsity and its final rank found 
% vectors, one for each matrix that is sliced 
minmse_miss = zeros(size(concentrations,2),1); 
minwmse_miss = zeros(size(concentrations,2),1);
minmse_fill = zeros(size(concentrations,2),1);
minRAD = zeros(length(concentrations),1);
minwmse_fill = zeros(size(concentrations,2),1);
min_fn = zeros(size(concentrations,2),1);
min_fn_fill = zeros(size(concentrations,2),1);
X_pred_best = zeros(dim1,dim2,size(concentrations,2)); % and the X predictions
R2 = zeros(size(concentrations,2),1);
% initialise error vars 
% column for each matrix (fn vs sparsity)
mse_miss = zeros(size(fns,2),size(concentrations,2));
smse_miss = zeros(size(fns,2),size(concentrations,2));
msefill= zeros(size(concentrations,2),1);
RADfill = zeros(size(fns,2),size(concentrations,2));
wmse_miss = zeros(size(fns,2),size(concentrations,2));
wmsefill = zeros(size(fns,2),size(concentrations,2));
RAD = zeros(size(fns,2),size(concentrations,2));
cumulative_contribution = zeros(length(concentrations),length(fns));
SVs = zeros(length(fns),length(concentrations));

%For gavish 
ranks_attempted = zeros(1000,length(concentrations));
%time the method 
tic
j=0; % counter for intervals
plotcount = 0;
for c = concentrations
     
    j=j+1; % for concentrations
    disp('Concentration of component1')
    disp(c)
    %Import sparse matrix for toy problem 
    Xs = X(:,:,j);
    Xfilled = Xs;
    dim = size(Xs);
    dim1 = dim(1);
    dim2 = dim(2);
    missing_ind = find(isnan(Xs));
    [row,col] = find(~isnan(Xs));
    filled_ind = find(~isnan(Xs));
    
    % Find the optimal rank for this matrix
    
    if Gavish==1
        [min_fn(j),minmse_fill(j),minwmse_fill(j),minRAD(j), X_pred_best,~,~,iter(j)] = solveGavish(Xs, dim1, dim2,1e-10,2, mixtures, c);
        %[fn,msefill,wmsefill,R2, X_pred,cutoff,SVs,iterations] = solveGavish(Xs, dim1, dim2, conv, iter)
        ranks_attempted(:,j) = min_fn(j);
    else
        % Iterative PCA with wMSE or MSE used to find rank

        %Find rank by minimising the mse or wmse 
        i=0;
        for fn=fns
            i=i+1;
            disp('rank')
            disp(fn)
            %[U,D,V,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
            [U,D,V,St,X_pred, iters]=missing_svd(Xs,fn,1,1,1e-10,100000,1);
            %X_pred is the model predictions, not only missing values are
            %filled 
            Xm = St*V';% the actual model, not just filled values - can only use errors of filled values to find model error
            Xm= X_pred;
            msefill(i,j) = (sum((Xm(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
            wmsefill(i,j) = find_wmse(Xs(filled_ind), Xm(filled_ind), length(filled_ind));
            abserrorfill = abs(X_pred(filled_ind)-Xs(filled_ind));
            RAD = abs(abserrorfill./Xs(filled_ind));
            RAD = RAD(find(~isnan(RAD)));
            RAD = RAD(find(~isinf(RAD)));
            RADfill(i,j) =sum(RAD)/length(RAD);
            %Cyt = corrcoef(X_pred(missing_ind),Xtrue(missing_ind));
            
            %R(i,j)=sqrt(Cyt(2,1));
        end
        
        minmse_fill(j)= min(msefill(:,j));
        min_fn_fill(j) = fns(find(msefill(:,j)==minmse_fill(j)));
        minwmse_fill(j) = min(wmsefill(:,j));
        minmse_miss(j) = min(mse_miss(:,j));
        minwmse_miss(j)=min(wmse_miss(:,j));
        
        if winsorized_mse ==1
            min_index = find(wmsefill(:,j)==minwmse_fill(j));
        else
            min_index = find(msefill(:,j)==minmse_fill(j));
        end 
        min_fn(j) = fns(min_index);
        %R2(j) = R(min_index,j)^2;
        minRAD(j) = RADfill(min_index,j);
        [U,D,V,St,X_pred_best(:,:,j), ~]=missing_svd(Xs,min_fn(j),1,1,1e-10,100000,1);
        % plot the scree plot for each composition
%         plotcount =plotcount+1;
%         subplot(6,2,plotcount)
%         [U,D,V,St,X_pred_plot, ~]=missing_svd(Xs,20,1,1,1e-3,1000);
%         y = diag(D);
%         SVs(:,j)=y;
%         semilogy(1:20, y(1:20))
%         hold on 
%         semilogy(min_fn(j), y(min_fn(j)),'ro')
%         xlabel('Rank')
%         ylabel('Singular values')
%         plotcount =plotcount+1;
%         subplot(6,2,plotcount)
%         plot(1:20,cumsum(diag(D))/sum(diag(D)))
%         cumulative_contribution(j,1:20) = cumsum(diag(D))/sum(diag(D));
%         xlabel('Rank')
%         ylabel('Cumulative contribution to error')
    end 
end 
% Create table with results 
Results = table(["Concentration of compound 1 (%)"; (concentrations'*100)],[ "Rank";min_fn],[ "SMSE (J/mol)";sqrt(minmse_fill)],[ "wSMSE (J/mol)";sqrt(minwmse_fill)], ["AARD (%)";minRAD*100]);
disp(Results)
filename = strcat('2wayRunPCASmall', date, '.mat');
save(filename)
toc


%% Error bars
error_boot = zeros(length(filled_ind),1);
for j=1:length(filled_ind)
    X_pred_best_boot(LOOCV_removed_row(j,c), LOOCV_removed_col(j,c)) = mean(X_pred_LOOCV(LOOCV_removed_row(j,c), LOOCV_removed_col(j,c), :));
    X_var_best_boot(LOOCV_removed_row(j,c), LOOCV_removed_col(j,c)) = var(X_pred_LOOCV(LOOCV_removed_row(j), LOOCV_removed_col(j,c), :));
    % actual error of this averaged prediction
    error_boot(j)= X_pred_best_boot(LOOCV_removed_row(j,c), LOOCV_removed_col(j,c))-Xs(LOOCV_removed_row(j,c), LOOCV_removed_col(j,c));
end
mse_boot = sum(error_boot.^2)/length(error_boot);
%% Plots of the singular values 
clf
dim2=40;
c = 0.5200;
if import ==1
    % import the relevant sparse matrix from the spreadsheet
    Ts = readtable(filename, 'Sheet', num2str(c));
    Xs = table2array(Ts);
    missing_ind = find(isnan(Xs));
    filled_ind = find(~isnan(Xs));
end 
ind = find(concentrations == c);

[U,D,V,St,X_pred_plot,~]=missing_svd(Xs,dim2,1,1,1e-5,1000); %uncomment to
mseplot = (sum((X_pred(missing_ind)-X(missing_ind)).^2))/length(remove_ind);
wmseplot= find_wmse(X(missing_ind), X_pred_plot(missing_ind), length(missing_ind));
smseplot = sqrt(mse_miss);
Cyt = corrcoef(X_pred(missing_ind),Xs(missing_ind));
R=sqrt(Cyt(2,1));

%choose rank
subplot(2,1,1)
y = diag(D);
semilogy(1:dim2, diag(D))
xlabel('Rank')
ylabel('Singular values')

if (rank_mat)>0
    hold on 
    semilogy(rank_mat, y(rank_mat), 'ko','MarkerSize',5, 'LineWidth',5)
    legend('SVs', 'True rank')
    hold off
end 
subplot(2,1,2)
plot(1:dim2,cumsum(diag(D))/sum(diag(D)))
xlabel('Rank')
ylabel('Cumulative contribution to error')
if rank_mat>0
    hold on
    y2 = y(1:rank_mat);
    realranky = cumsum(y(1:dim2));
    plot(min_fn(ind), realranky(min_fn(ind)-fns(1)+1)/sum(diag(D)), 'r*', 'MarkerSize',15, 'LineWidth',1.5)
    legend('Cumulative error', 'Predicted rank')
    hold off
end 
sgtitle(strcat('Sparsity = ',num2str(c)))
%% Plots of errors 
clf
plotno = length(concentrations);
a = ceil(sqrt(plotno));
b = ceil(plotno/a);
i=0;
for c = concentrations
    i=i+1;
    ind = find(concentrations == c);
    subplot(a,b,i)
    titlestr = strcat('Sparsity = ',num2str(c));
    plot(fns, wmse_miss(:,ind))
    hold on
    plot(fns,mse_miss(:,ind), 'm')
    plot(rank_mat,mse_miss(rank_mat-fns(1)+1,ind), 'ko', 'MarkerSize',5, 'LineWidth',5)
    plot(min_fn(ind),wmse_miss(min_fn(ind)-fns(1)+1,ind), 'r*', 'MarkerSize',15, 'LineWidth',1.5)
    hold off
    legend('Winsorized MSE', 'MSE', 'True rank', 'Predicted Rank')
    title(titlestr)
    xlabel('Rank')
    ylabel('Error')
end 
%% Plot for Gavish 
if Gavish ==1 
    clf
    subplot(2,1,1)
    plot(1:dim2,d_original)
    xlabel('Rank')
    ylabel('Singular value')
    xlabel('Component number')
    hold on 
    plot(rank_mat, d_original(rank_mat), 'ko')
    plot(fn, d_original(fn), 'ro')
    legend('Singular values', 'True rank', 'Predicted rank')
    
    subplot(2,1,2)
    plot(1:dim2,cumsum(d_original)/sum(d_original))
    xlabel('Component number')
    ylabel('Cumulative contribution to error')
end 

%% Functions

function  [fn,msefill,wmsefill,RAD, X_pred,cutoff,SVs,iterations] = solveGavish(Xs, dim1, dim2, conv, iter, mixtures, conc)
%[min_fn(j),minmse_fill(j),minwmse_fill(j),R2(j), X_pred_best,ranks] = solveGavish(Xs, dim1, dim2,1e-3,1000)
    % Gavish Hard thresholding used to find the rank of the matrix using PCA
    % with SVD 
    % Input 
    % Xs = data matrix 
    % dim1 = number of rows
    % dim2 = number of columns 
    % conv = convergence criterion 
    % iter = maximum number of iterations 
    % Output 
    % fn = Number of factors/ PCs
    % mse = mean squared error for observed entries compared to the model predictions for these for the optimal number of factors 
    % wmse = winsorized mean squared error for observed entries compared to the model predictions for these for the optimal number of factors
    % R2 = correlation between the observed entries and the model predictions
    % of those entries
    % X_pred = the data matrix as predicted by the optimal number of factors 
    % cutoff = cutoff used to find the number of factors 
    % SVs = singular value matrix found 
    filled_ind = find(~isnan(Xs));
    % find omega 
    if dim1/dim2 ==1
        omega = 2.858; % omega(beta)=2.858 for n*n square matrix
    elseif dim1/dim2 < 1
        beta = dim1/dim2;
        omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
    else
        beta = dim2/dim1; 
        omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
    end 
    % PCA done once for the maximum number of factors, the matrix needs to be filled to do this, which is
    % done in the missing_svd function 
    %Xfilled = fill_data(Xs);
    [Xfilled, missing] = filldata3(Xs, 'uni',mixtures,conc, 'none', 298.15);
    if dim1<dim2
        fn = dim1;
    else 
        fn=dim2;
    end 

    % find SVs for matrix 
    [U,SVs,D,~,~,~, ~]=missing_svd(Xfilled,fn,1,1,1e-4,2,1);
    %[U,D,V,St,X_pred, AllSVs, iterations]=missing_svd(X,fn,center,scale,conv,max_iter, use_missing)
    % find the threshold 
    y_med = median(diag(SVs)); %ymed = median singular value of the noisy matrix  
    cutoff = omega*y_med; %cutoff= tau= omega(beta)*ymed; matrix
     % Keep modes w/ sig > cutoff; rank chosen as hard cutoff
    fn = length(find(diag(SVs)>cutoff));
    % reconstruct with the new SVs
    SVs(:,fn:end)=0;
        
    % [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd(X,fn,center,scale,conv,max_iter)
     % solve for the correct number of factors 
    [~,~,D,St,X_pred,~,iterations] = missing_svd(Xs,fn,1,1,conv,iter,1);
    msefill = (sum((X_pred(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
    wmsefill = find_wmse(Xs(filled_ind), X_pred(filled_ind), length(filled_ind));
    RAD = (abs((X_pred(filled_ind)-Xs(filled_ind))./Xs(filled_ind)));
    RAD = RAD(find(~isnan(RAD)));
    RAD = RAD(find(~isinf(RAD)));
    RAD = sum(RAD)/length(RAD);
    Cyt = corrcoef(X_pred(filled_ind),Xs(filled_ind));
    R2=(Cyt(2,1));
end 

function [X_filled, missing] = filldata2(X, method,mixtures,concintervalarray)
    % filldata3 fills a 3-way array with the arithmetic mean of the other
    % values in the array
    % Input 
    % X = data array 
    % method = 'avg' - column and row averages 
    %       or 'uni' - unifac predictions 
    % components = components ordered according to their place on the axis 
    % concinterval = concentration interval over which to fill data 
    % Output 
    % X_filled = filled data array
    % missing = linear indices of missing data 

    % fill as done for PCA with SVD (averages of each dimension)
    [i, j]=findnan3(X);% returns rows and columns with nan elements
    dim = size(X);
    missing = isnan(X); %linear indices 

    X_filled = X;
    % fill missing values with zeros 
   
    
    if strcmp(method,'avg')
    % find all means 
        mean1 = sum(X_filled,1)./(ones(1,dim(2))*dim(1)-sum(missing,1)); % (1,j,k) dim1 
        mean2= sum(X_filled,2)./(ones(dim(1),1)*dim(2)-sum(missing,2)); % (i,1,k)  dim2 
        %replace nan means with 0 
        mean1(find(isnan(mean1)))=0;
        mean2(find(isnan(mean2)))=0;
        for ind =1:length(i)
           % for all NaN elements that exist, loop through them to replace
            X_filled(i(ind),j(ind))=(mean1(1,j(ind))+mean2(i(ind),1))/2;
        end      
    else %method ='uni'
        % linearise X to a column vector - entries match mixtures
        X1 = reshape(X(:,:,1),[dim(1)*dim(2),1]); % column vector used for checks and indices
        Xtemp = reshape(X,[dim(1)*dim(2),1]); % reshaped to column vectors with each column containing a concentration interval
        
        mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
        mix1 = mixtures(:,[1,2]);
        mix2 = mixtures(:,[3,4]);
        % load prediction data and mixtures  
        load('heUNIQUACforT=298.15.mat','he') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
        load('heUNIQUACforT=298.15.mat','mixture') % mixtures for he data - 5151 components
        load('heUNIQUACforT=298.15.mat','conc_interval')
        %convert arrays of concentrations to strings to allow them to be
        %compared 
        conc_unifac = string(conc_interval);
        conc_array = string(concintervalarray);
        conc_array2 = string(1-conc_interval);
        for c = 1:length(conc_array)
            concindices(c) = find(strcmp(conc_array(c),conc_unifac));
            concindices2(c) = find(strcmp(conc_array2(c),conc_unifac));
        end 
        %components = possible components in this array, ordered 
        for ind = 1:length(X1)
            if isnan(X1(ind,1))
                % fill with UNIFAC prediction 
                [~,indexpred] = ismember(mixturesarray(ind,:),mixture,'rows');
                
                if indexpred ==0
                    %mixture could be swapped around ie the concentration
                    %is 1-conc as well and is not in the UNIFAC data set 
                    [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixtures,'rows');
                    if indexpred == 0
                        Xtemp(ind,:) = 0;
                        disp(ind)
                    else 
                        Xtemp(ind,:)=he(indexpred,concindices2);
                    end    
                else
                %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                    Xtemp(ind,:)=he(indexpred,concindices);
                end 
            end 
        end 
       
        %reshape X to the 3-way array
        X_filled = reshape(Xtemp,[dim(1), dim(2)]);
        disp(X_filled)
        %fills remaining few Nan values with averages 
        X_filled = filldata2(X_filled,'avg',mixtures,conc_interval);
    end 
end 

function [i,j]=findnan3(X)
% Input 
% X = data array 
% Output 
% i, j, k = indexes of nan values in X

% isnan(X) returns logical array 
% all(isnan(X),1) is also a logical array (1x30x40) - each row
% find(all(isnan(X),1)) returns linear indices -> reshape array to correct
% dimension to find missing slabs 
% findnan3 finds the nan indices of a 3-way array 
    

    i=[];
    j=[];
    
    
        % done per z slice 
        [itemp, jtemp]= find(isnan(X));
        i = [i; itemp];
        j = [j;jtemp];
    
end

function [indices]=reorder_He(X, rand)
% X = matrix to be ordered 
% indices = new indices to be used 
% rand = manner of reordering 
    %rand = 'any' = reorder in any way, default 
    %rand = 'rows' = reorder values within rows 
    % rand = 'cols' = reorder values within columns 
    % rand = 'zeroslineone' = reorder values so that the diagonal ends up
    % on the 1st row, not done yet 
    
    if rand =="any"
        indices = randperm(numel(X));
    elseif rand == "rows" 
        indices = 1:1:numel(X);
        n = size(X,1); %no rows
        m = size(X,2); %no cols
        for i =1:n
            p = (i-1)*m+1:(i)*m;
            temp = indices(p);
            temp = temp(randperm(length(temp)));
            indices(p) = temp;
        end 
    elseif rand == "cols"
        indices = 1:1:numel(X);
        n = size(X,1); %no rows
        m = size(X,2); %no cols
        for i =1:m
            p = (i-1)*n+1:n:(n*m)-(n-1);
            temp = indices(p);
            temp = temp(randperm(length(temp)));
            indices(p) = temp;
        end
        
    else 
        indices = randperm(numel(X));
    end 
end 

function [X,Xnoise, rankU, rankV, noise]=create_matrix(n,m,r, mu, sigma)
    % create a matrix of size nxm with rank =r using SVD 
    % Input 
    % n = dimension 1 / rows
    % m = dimension 2 / columns
    % r = rank of matrix 
    % mu = mean of data 
    % sigma = standard deviation of data 
    % Output 
    % X = data matrix formed; X = UDV'
    % Xnoise = X+noise; noisy matrix  
    % rankU = rank of U 
    % rankV = rank of V
    % noise = noise matrix, normally distributed 
    
    sing_vals = normrnd(mu, sigma, r,1); % generate r singular values that are normally distributed 
    sing_vals = sort(sing_vals, 'descend'); % sorted from highest to lowest
    D = zeros(n,m);%initialise D 
    for i = 1:r
        D(i,i) = sing_vals(i); %populate array with values 
    end 
    % fill the matrices U,V with random entries
    U = randn(n,n);
    V = randn(m,m);
    % find orthonormal bases of U and V 
    rankU = rank(U);
    rankV = rank(V);
    U = orth(U);
    V = orth(V);
    %fill the final matrix 
    X = U*D*V'; %create the matrix of rank r
    noise = randn(size(X))/(n*m);
    Xnoise = X+noise;
end 

function [X_sparse,M, remove_ind, fill_ind]=fill_matrix(X, perc_remove)
    % Fill data into a matrix until you have reached 1-perc_remove
    % Input 
    % X = filled matrix 
    % perc_remove = percentage of the data to remove. Pass as a fraction 
    % Output 
    % X_sparse = sparse matrix 
    % M = matrix with the same dimensions of X with a 1 for
    % each observed entry in X_sparse and a 0 for an unobserved entry
    % remove_ind = linear indices of unobserved entries 
    % fill_ind = linear indices of observed entries 
    
    [n,m]=size(X);
    num_removed = floor(m*n*perc_remove);
    %create NaN matrix 
    X_sparse = NaN(n,m);
    M = zeros(n,m);
    
    % fill in one entry in every row and column, randomly 
    filled_rows = zeros(n,1);
    cols=randperm(m); % sort the columns in a random order 
    for j =cols
        i = randi([1,n],1);
        while ismember(i,filled_rows)
            i = randi([1,n],1);
        end
        filled_rows(j)=i;
        X_sparse(i, j)=X(i,j);
        M(i,j)=1;
    end 
    
    for k= 1:(m*n-num_removed-m) % add this many values from array X, 
        i = randi(n,1); 
        j = randi(m,1);
        while M(i,j)==1
            %repeat the finding of indices to fill if this index is already filled 
            i = randi(n,1); %row
            j = randi(m,1);%col
        end 
        X_sparse(i,j)= X(i,j); %reassign value as original X value 
    end 
    remove_ind = find(isnan(X_sparse));
    fill_ind = find(~isnan(X_sparse));
end 

function [X_sparse,remove_ind,fill_ind]=remove_matrix(X,perc_remove)
% Remove data from the X matrix until the percentage is removed, must
% remove entire rows or columns for the matrix to make sense
    [m,n]=size(X);
    num_removed = floor(m*n*perc_remove);
    
    X_sparse = X;
    for k= 1:(num_removed) % remove this many values from array X 
        i = randi(m,1); %matrix inidices start at 1 in MATLAB
        j = randi(n,1);
        while isnan(X_sparse(i,j))
            i = randi(m,1); %repeat the removal if necessary  
            j = randi(n,1);
        end 
        X_sparse(i,j)= nan; %reassign value as NaN value 
    end 
    remove_ind = find(isnan(X_sparse));
    fill_ind = find(~isnan(X_sparse));
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

function [X_filled]=fill_data(X)
% Fill a matrix which has missing entries. Missing entries are each filled with the average of the observed entries in its row and column. 
% Input 
% X = matrix with missing data 
% Output 
% X = filled matrix 
    [m,n]=size(X);
    missing_ind = (isnan(X));
    [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
    X_filled=X;
    X_filled(missing_ind)=0; %fill NaN values with 0
    mean_col = sum(X_filled,1)./(ones(1,n)*m-sum(missing_ind,1)); 
    mean_row = sum(X_filled,2)./(ones(1,m)*n-sum(missing_ind,2));  
    % for all NaN elements that exist, loop through them to replace with means 
    for k =1:length(i) 
        X_filled(i(k),j(k))=(mean_row(i(k))+mean_col(j(k)))/2;
    end  
end 

function [U,D,V,St,X_pred, AllSVs, iterations]=missing_svd(X,fn,center,scale,conv,max_iter, use_missing,methodguess)
    % Fill a matrix of missing data using PCA with SVD and a given number of
    % PCs. Can also handle non-missing data. Missing data is handled as NaN
    % values 
    %
    % Input 
    % X = matrix
    % fn = number of factors/ PCs used to fill the matrix 
    % center = 1 center data / = 0 do not center 
    % scale = 1 scale data / = 0 do not scale
    % conv = stopping criterion, absolute value of the relative change of the
    % sum of squares of the values of the unobserved entries 
    % max_iter = maximum number of iterations 
    % use_missing = use the missing entries for the convergence, =1 to use
    % missing 
    % Output 
    % S,V,D from X = SVD'
    % St = SV
    % X_pred = filled X with new values for the missing entries
    [m,n]=size(X);
    missing_ind = find(isnan(X));
    filled_ind = find(~isnan(X));
    %choose indices to use for convergence 
    if use_missing ==1
        indices=missing_ind;
    else 
        indices = filled_ind;
    end 
     
    if any(isnan(X)) % there is missing data 
        Xfilled = fill_data(X); 
        SS = sum(sum(Xfilled(indices).^2));
        f=2*conv;
        iter = 0;
        while iter<max_iter && f>conv
            iter = iter+1;
            SSold = SS; 
            % preprocess = scale and then center 
            mx = mean(Xfilled);
             
            if center ==1
                mx = mean(Xfilled,2); %columnwise mean 
                Xc = Xfilled-ones(1,n)*mx; %centering of the data done - subtract mean from each entry in a column 
            else 
                Xc=Xfilled;
            end 
            if scale ==1 
                sj = sqrt(sum(sum((Xc).^2)));
                Xc = Xc/sj;
            end
            %perform SVD
            [U,D,V]=svd(Xc);
            AllSVs = D;
            St=U*D;
            St = St(:,1:fn);
            D=D(:,1:fn); %same as D(:,fn:end)=0 and others kept the same for multiplication
            %D(:,fn:end)=0;
            U=U(:,1:fn);
            V=V(:,1:fn);
            % post process = uncenter, unscale  
            X_pred = St*V';
            if scale ==1
                X_pred = X_pred*sj;
            end
            if center ==1
                X_pred =X_pred+ones(1,n)*mx;
            end 
            % fill misssing values 
            Xfilled(missing_ind)=X_pred(missing_ind);
            SS = sum(sum((X_pred(indices)).^2));
            %f = (sum(X_pred(filled_ind)-X(filled_ind).^2)); %residuals -
            %does not work 
            f = abs(SS-SSold)/(SSold);
        end
        iterations = iter;
        
    else 
        iterations =1;
        % no missing data 
        Xfilled=X;
        if scale ==1 
            sj = sqrt(sum(sum((Xfilled).^2)));
            Xfilled = Xfilled/sj;
        end
        if center ==1
            mx = mean(Xfilled);
            %Xc = normalize(Xf); %normalizing 
            Xc = Xfilled-ones(m,1)*mx; %centering
        else 
            Xc=X;
        end 
        [U,D,V]=svd(Xc);
        St=U*D;
        AllSVs=D;
        St = St(:,1:fn);
        D=D(:,1:fn);
        U=U(:,1:fn);
        V=V(:,1:fn);
         
        if center ==1
            X_pred = St*V'+ones(m,1)*mx;
        else 
            X_pred = St*V';
        end 
        if scale ==1
                X_pred = X_pred*sj;
        end
    end %end if  else 
end % end function 

function [wmse]=find_wmse(X,X_pred,no_points)
% Find the winsorized mse between the matrices X and X_pred 
% Input 
% X = true matrix 
% X_pred = predicted matrix using the model 
% both X and X_pred must be passed as vectors 
% no_points = number of points in the vectors passed
% Output 
% wmse = the winsorized MSE

    perc5 = prctile(X_pred,5, 'all');
    perc95 = prctile(X_pred, 95, 'all');
    %reassign
    X_pred(X_pred<perc5)=perc5;
    X_pred(X_pred> perc95)=perc95;
    wmse = (sum((X_pred-X).^2))/no_points;
end 

function [wmse]=find_wmse_error(errors, count)
    perc5=prctile(errors,5,'all');
    perc95=prctile(errors,95,'all');
    errors(errors<perc5)=perc5;
    errors(errors>perc95)=perc95;
    wmse = (sum((errors).^2))/count;
end 
