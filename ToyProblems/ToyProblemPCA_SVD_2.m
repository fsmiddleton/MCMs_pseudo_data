%% Toy problem to create code for SVD used for matrix completion 
% Francesca Middleton, 2022-03-02

%% Import or create data matrix for testing 
clc
clear

import = 1; 
% import = 0 to create a low rank matrix, 
% import = 1 to fetch an already created low rank matrix (2-way) 
% import = 2 to import excess enthalpy data (2-way),

filename = 'RandomMatrixNoNoise.xlsx';
export = 0; % either export data(1) or don't, only used for creating a low rank matrix
noise=1; %noise =1 to add noise to the creation of the matrix
rank_mat = 0; % to be assigned later 

if import ==0
    % create array
    % specify size and rank of array, choosing random mu and sigma to create
    % singular values from, and the noise  
    dim1=50;
    dim2=50;
    rank_mat = 10;
    mu = 20;
    sigma = 3;
    [Xs,Xnoise, rankU, rankV, noiseMat]=create_matrix(dim1,dim2,rank_mat, mu, sigma);
    Xdiff=Xs-Xnoise;
    if noise ==1
        X=Xnoise;
    else 
        X = Xs;
    end 
     
    if export ==1
        %create table with all this information
        Ts = array2table(Xs);
        Tnoise = array2table(Xnoise);
        Tfinal = array2table(X);
        info = table(["rank_mat"; "mu"; "sigma"; "noise"],[rank_mat; mu; sigma; noise], ["Sheet1";"Sheet2";"Sheet3";"Sheet4"], ["Info"; "Random matrix";"Noise matrix";"Noisy matrix"] );
        %export it 
        sheetname = append('Rank', num2str(rank_mat));
        writetable(info,filename,'Sheet',1)
        writetable(Ts,filename,'Sheet',2)
        writetable(Tnoise,filename,'Sheet',3)
        writetable(Tfinal,filename,'Sheet',4)
        for i = 0.5:0.01:0.6
            [Xsparse,missing_ind,filled_ind]=fill_matrix(X,i);
            Tsparse = array2table(Xsparse);
            writetable(Tsparse, filename, 'Sheet', num2str(i))
        end 
    end
    
        
elseif import ==1
    %import the data from the previously created random matrix 
    filename = 'RandomMatrixNoise.xlsx'; % choose file here
    % matrix with and without noise in sheets 2 and 4
    if noise == 1
        T2 = readtable(filename, 'Sheet', 4);
    else
        T2 = readtable(filename, 'Sheet', 2);
    end 
    X = table2array(T2);
    T3 = readtable(filename, 'Sheet', 1); %info in sheet 1
    rank_mat = T3(1,2);
    
elseif import ==2
    % Import excess enthalpy data for one matrix 
    T1 = readtable('HEmatrix298.15.2.xlsx', 'Sheet', 3); 
    X = table2array(T1);%only use T1, small enough to manage
    
else
    disp('Please choose 0, 1, 2 to use data')
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
%% Find optimal rank of the model for the filled data 

clc
clear

%Import true data 
filename = 'RandomMatrixNoise.xlsx';
Ttrue = readtable(filename, 'Sheet', 4);
Xtrue = table2array(Ttrue);
% Import other data for most missing to allow variables to be declared
miss=30;
T = readtable(filename, 'Sheet', num2str(miss/100));
Xs=table2array(T);

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
winsorized_mse = 1; %1=use wmse. Used for iterative PCA 
 % Choose sparsities used for finding the rank 
missing = 30:10:80;
%ranks to try 
fns = 1:1:20;
%choose whether to reorder matrix or not, 1 = true 
reorder = 0;
        
%intialise the metrics to analyse each sparsity and its final rank found 
% vectors, one for each matrix that is sliced 
minmse_miss = zeros(size(missing,2),1); 
minwmse_miss = zeros(size(missing,2),1);
minmse_fill = zeros(size(missing,2),1);
minwmse_fill = zeros(size(missing,2),1);
min_fn = zeros(size(missing,2),1);
min_fn_fill = zeros(size(missing,2),1);
X_pred_best = zeros(dim1,dim2,size(missing,2)); % and the X predictions
R2 = zeros(size(missing,2),1);
% initialise error vars 
% column for each matrix (fn vs sparsity)
mse_miss = zeros(size(fns,2),size(missing,2));
smse_miss = zeros(size(fns,2),size(missing,2));
msefill= zeros(size(missing,2),1);
wmse_miss = zeros(size(fns,2),size(missing,2));
wmsefill = zeros(size(fns,2),size(missing,2));
R = zeros(size(fns,2),size(missing,2));
cumulative_contribution = zeros(length(missing),length(fns));
SVs = zeros(length(fns),length(missing));

%time the method 
tic
j=0; % counter for intervals
plotcount = 0;
for miss = missing
     
    j=j+1; % for sparsity
    disp('missing')
    disp(miss)
    %Import sparse matrix for toy problem 
    T = readtable(filename, 'Sheet', num2str(miss/100));
    Xs = table2array(T);
    Xs = remove_nan2(Xs);
    Xfilled = Xs;
    dim = size(Xs);
    dim1 = dim(1);
    dim2 = dim(2);
    missing_ind = find(isnan(Xs));
    [row,col] = find(~isnan(Xs));
    filled_ind = find(~isnan(Xs));
    
    % Find the optimal rank for this matrix
    
    if Gavish==1
        [min_fn(j),minmse_miss(j),minwmse_miss(j),R2(j), X_pred_best(:,:,j)] = solveGavish(Xs, dim1, dim2,1e-3,200);
         
    else
        % Iterative PCA with wMSE or MSE used to find rank

        %Find rank by minimising the mse or wmse 
        i=0;
        for fn=fns
            i=i+1;
            disp('rank')
            disp(fn)
            %[U,D,V,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
            [U,D,V,St,X_pred, ~]=missing_svd(Xs,fn,1,1,1e-3,1000);
            
            SRSS = sqrt((sum((X_pred(missing_ind)-Xtrue(missing_ind)).^2)));
            mse_miss(i,j) = (sum((X_pred(missing_ind)-Xtrue(missing_ind)).^2))/length(missing_ind);
            wmse_miss(i,j)= find_wmse(Xtrue(missing_ind), X_pred(missing_ind), length(missing_ind));
            smse_miss(i,j) = sqrt(mse_miss(i));
            Xm = St*V';% the actual model, not just filled values - can only use errors of filled values to find model error
            Xm= X_pred;
            msefill(i,j) = (sum((Xm(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
            wmsefill(i,j) = find_wmse(Xs(filled_ind), Xm(filled_ind), length(filled_ind));
            Cyt = corrcoef(X_pred(missing_ind),Xtrue(missing_ind));
            R(i,j)=sqrt(Cyt(2,1));
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
        R2(j) = R(min_index,j)^2;
        [U,D,V,St,X_pred_best(:,:,j), ~]=missing_svd(Xs,min_fn(j),1,1,1e-3,1000);
        % plot the scree plot for each composition
        plotcount =plotcount+1;
        subplot(6,2,plotcount)
        [U,D,V,St,X_pred_plot, ~]=missing_svd(Xs,20,1,1,1e-3,1000);
        y = diag(D);
        SVs(:,j)=y;
        semilogy(1:20, y(1:20))
        hold on 
        semilogy(min_fn(j), y(min_fn(j)),'ro')
        xlabel('Rank')
        ylabel('Singular values')
        plotcount =plotcount+1;
        subplot(6,2,plotcount)
        plot(1:20,cumsum(diag(D))/sum(diag(D)))
        cumulative_contribution(j,1:20) = cumsum(diag(D))/sum(diag(D));
        xlabel('Rank')
        ylabel('Cumulative contribution to error')
    end 
end 
% Create table with results 
Results = table(["Sparsity"; (missing')],[ "SVs";min_fn],[ "MSE";minmse_fill],[ "wMSE";minwmse_fill], ["MSE missing";minmse_miss]);
disp(Results)
toc
%% LOOCV to find rank using mse  
clc
clear


%Import true data 
filename = 'RandomMatrixNoise.xlsx';
Ttrue = readtable(filename, 'Sheet', 4);
Xtrue = table2array(Ttrue);
% Import other data for most missing to allow variables to be declared
miss=30;
T = readtable(filename, 'Sheet', num2str(miss/100));
Xs=table2array(T);


%[X_3D, Xtrue_3D, Factors] = importX3D(filename, dim);
%unfolding

%Xs = reshape(X_3D,[],40); % columns remain for unfolding 

%time the code 
tic
% declare missing % used and the file from which to extract information 
missing = 30:10:80;
%ranks to try 
fns = 1:1:20; % maximum is 50 
%necessary vars 
dim = size(Xs);
dim1=dim(1);
dim2=dim(2);
missing_ind = find(isnan(Xs));
% largest length of filled_linear_ind is for 0.3
filled_ind = find(~isnan(Xs));

% choose whether to use mse or wmse to find the optimal rank 
winsorized_mse = 0;
% declare vars for analysis 
mse_LOOCV = zeros(length(missing),length(fns));
wmse_LOOCV = zeros(length(missing),length(fns));
mse_LOOCV_missing = zeros(length(missing),length(fns));
min_mse_LOOCV = zeros(length(missing),1);
min_wmse_LOOCV = zeros(length(missing),1);
min_mse_LOOCV_missing = zeros(length(missing),1);
X_pred_LOOCV = zeros(dim1,dim2,length(filled_ind));
LOOCV_removed_col = zeros(length(filled_ind),1);
LOOCV_removed_row = zeros(length(filled_ind),1);
SVs = zeros(dim1,length(missing));

% Error bars - LOOCV 
% import the relevant sparse matrix from the spreadsheet
count_missing = 0;
for miss = missing
    count_missing = count_missing+1;
    T = readtable(filename, 'Sheet', num2str(miss/100));
    Xs = table2array(T);
    Xs = remove_nan2(Xs);
    Xfilled = Xs;
    dim = size(Xs);
    dim1 = dim(1);
    dim2 = dim(2);
    missing_ind = find(isnan(Xs));
    [row,col] = find(~isnan(Xs));
    filled_ind = find(~isnan(Xs));
    disp('% missing')
    disp(miss)
    %declare vars with size dependent on the array used 
    Xm_boot=zeros(length(filled_ind),1);
    error_LOOCV = zeros(length(filled_ind),1);
    mse_miss = zeros(length(missing_ind),1);
    %loop through ranks
    disp('fn')
    for fn = fns 
        disp(fn)
        k=0;% for the bootstrapping/ LOOCV for each 
        for filled_ind = filled_ind' %filled_linear_ind must be a row vector for a for loop    
            k=k+1;
            % remove a point from Xs
            X_b = Xs;
            X_b(filled_ind) = nan;
            if find(~isnan(X_b(:,col(k)))) & find(~isnan(X_b(row(k),:))) %ensure at least one value in each column and row
                LOOCV_removed_col(k,count_missing) = col(k);
                LOOCV_removed_row(k,count_missing) = row(k);
                %perform iterative PCA on the slightly more empty matrix 
                %  [S,V,D,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
                [U,D,V,St,Xpred, ~]=missing_svd(X_b,fn,1,1,1e-3,1000);
                
            end
            Xm= St*V';
            error_LOOCV(k) = Xs(filled_ind)-Xm(filled_ind);
            mse_miss(k) = sum((Xtrue(missing_ind)-Xm(missing_ind)).^2)/length(missing_ind);
            Xm_boot(k) = Xm(filled_ind);
        end
        % mse for this composition and rank 
        mse_LOOCV(count_missing,fn)= sum(error_LOOCV.^2)/length(error_LOOCV);
        wmse_LOOCV(count_missing, fn) = find_wmse(Xs(filled_ind), Xm_boot, length(filled_ind));
        % consider how it did for the other missing values 
        % average of mses
        mse_LOOCV_missing(count_missing, fn) = sum(mse_miss)/length(mse_miss);
    end
    
    % find the optimal rank 
    if winsorized_mse ==1
        fn_LOOCV(count_missing) = find(wmse_LOOCV(count_missing,:) == min(wmse_LOOCV(count_missing,:)));
    else 
        fn_LOOCV(count_missing) = find(mse_LOOCV(count_missing,:) == min(mse_LOOCV(count_missing,:)));
    end
    [U,D,V,St,Xpred, SVsmiss]=missing_svd(Xs,fn_LOOCV(count_missing),1,1,1e-3,1000);
    y = diag(SVsmiss);
    SVs(:,count_missing)=y;
    min_mse_LOOCV(count_missing) = mse_LOOCV(count_missing,fn_LOOCV(count_missing));
    min_wmse_LOOCV(count_missing) = wmse_LOOCV(count_missing,fn_LOOCV(count_missing));
    min_mse_LOOCV_missing(count_missing) = mse_LOOCV_missing(count_missing,fn_LOOCV(count_missing));
    % use the best rank to find error bars 
end 
Results = table(["Sparsity"; (missing')],[ "SVs";fn_LOOCV'],[ "MSE";min_mse_LOOCV],[ "wMSE";min_wmse_LOOCV], ["Missing data MSE";min_mse_LOOCV_missing]);
disp(Results)
toc 
%% Gavish to find rank  
clc
clear
% check rank found using Gavish 
z=1; % which matrix to use 
missing = 30:10:80;

%Import true data 
filename = 'RandomMatrixNoNoise.xlsx';
Ttrue = readtable(filename, 'Sheet', 4);
Xtrue = table2array(Ttrue);
% Import other data for most missing to allow variables to be declared
miss=30;
T = readtable(filename, 'Sheet', num2str(miss/100));
Xs=table2array(T);
Xs = remove_nan2(Xs);
dim = size(Xs);
dim1 = dim(1);
dim2 = dim(2);
fn_LOOCV = zeros(length(missing),1);
fn_Gavish = zeros(length(missing),1);
SVs = zeros(dim1,length(missing));
msemiss = zeros(length(missing),1);
msefill= zeros(length(missing),1);
wmsemiss = zeros(length(missing),1);
wmsefill = zeros(length(missing),1);
fnGavish = zeros(length(missing),1);

count_missing = 0;
for miss = missing 
    disp('%missing')
    disp(miss)
    count_missing = count_missing +1;
    % Import matrix 
    T = readtable(filename, 'Sheet', num2str(miss/100));
    Xs = table2array(T);
    Xs = remove_nan2(Xs);
    missing_ind = find(isnan(Xs));
    filled_ind = find(~isnan(Xs)); 
    % Find rank using Gavish 
    if dim1/dim2 ==1
        omega = 2.858; % omega(beta)=2.858 for n*n square matrix
    elseif dim1/dim2 < 1
        beta = dim1/dim2;
        omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
    else
        beta = dim2/dim1; 
        omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
    end 
    % finds all SVs using iterative PCA using smaller of dim1 and dim2 
    if dim1<dim2
        %only one iteration of SVD decomposition used 
        [U,D,V,St,Xpred_Gavish,SVs, iterations]=missing_svd(Xs,dim1,1,1,1e-3,1); 
    else 
        [U,D,V,St,Xpred_Gavish,SVs, iterations]=missing_svd(Xs,dim2,1,1,1e-3,1);
    end 
    % D = singular values 
    % find cutoff for these 
    % use cutoff to fill 
    % repeat until covergence 
    for i=1:50
        y_med = median(diag(D)); %ymed = median singular value of the noisy matrix  
        cutoff = omega*y_med; %cutoff= tau= omega(beta)*ymed; matrix
        % Keep modes w/ sig > cutoff; rank chosen as hard cutoff
        d=diag(D);
        fnGavish(count_missing) = length(find(d>cutoff));
    end 
    SVs(:,count_missing)= diag(D);
    d(1+fnGavish(count_missing):end)=0;
    [~,~,V,St,X_pred,~, ~] = missing_svd(Xs,fnGavish(count_missing),1,1,1e-3,1000);
    Xm = St*V';
    msefill(count_missing,1) = (sum((Xm(filled_ind)-Xtrue(filled_ind)).^2))/length(filled_ind);
    wmsefill(count_missing,1) = find_wmse(Xtrue(filled_ind), Xm(filled_ind), length(filled_ind));
    msemiss(count_missing,1) = (sum((X_pred(missing_ind)-Xtrue(missing_ind)).^2))/length(missing_ind);
    wmsemiss(count_missing,1) = find_wmse(Xtrue(missing_ind), X_pred(missing_ind), length(missing_ind));
end 
% end Gavish
%% Error bars
error_boot = zeros(length(filled_ind),1);
for j=1:k
    X_pred_best_boot(LOOCV_removed_row(j,count_missing), LOOCV_removed_col(j,count_missing)) = mean(X_pred_LOOCV(LOOCV_removed_row(j,count_missing), LOOCV_removed_col(j,count_missing), :));
    X_var_best_boot(LOOCV_removed_row(j,count_missing), LOOCV_removed_col(j,count_missing)) = var(X_pred_LOOCV(LOOCV_removed_row(j), LOOCV_removed_col(j,count_missing), :));
    % actual error of this averaged prediction
    error_boot(j)= X_pred_best_boot(LOOCV_removed_row(j,count_missing), LOOCV_removed_col(j,count_missing))-Xs(LOOCV_removed_row(j,count_missing), LOOCV_removed_col(j,count_missing));
end
mse_boot = sum(error_boot.^2)/length(error_boot);
%% Plots of the singular values 
clf
dim2=40;
miss = 0.5200;
if import ==1
    % import the relevant sparse matrix from the spreadsheet
    Ts = readtable(filename, 'Sheet', num2str(miss));
    Xs = table2array(Ts);
    missing_ind = find(isnan(Xs));
    filled_ind = find(~isnan(Xs));
end 
ind = find(missing == miss);

[U,D,V,St,X_pred_plot,~]=missing_svd(Xs,dim2,1,1,1e-3,1000); %uncomment to
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
sgtitle(strcat('Sparsity = ',num2str(miss)))
%% Plots of errors 
clf
plotno = length(missing);
a = ceil(sqrt(plotno));
b = ceil(plotno/a);
i=0;
for miss = missing
    i=i+1;
    ind = find(missing == miss);
    subplot(a,b,i)
    titlestr = strcat('Sparsity = ',num2str(miss));
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

function  [fn,mse,wmse,R2, X_pred,cutoff,SVs] = solveGavish(Xs, dim1, dim2, conv, iter)
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

    missing_ind = find(isnan(Xs));
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
    Xfilled = fill_data(Xs);
    SS = sum(sum(Xfilled(missing_ind).^2));
    f=conv*2;
    it=0;
    if dim1<dim2
        fn = dim1;
    else 
        fn=dim2;
    end 
    
    while f>conv && it<iter
        SSold = SS;
        % find SVs for matrix 
        [U,SVs,D,~,~,~, ~]=missing_svd(Xfilled,dim2,1,1,1e-4,2);
        % find the threshold 
        y_med = median(diag(SVs)); %ymed = median singular value of the noisy matrix  
        cutoff = omega*y_med; %cutoff= tau= omega(beta)*ymed; matrix
         % Keep modes w/ sig > cutoff; rank chosen as hard cutoff
        fn = length(find(diag(SVs)>cutoff));
        % reconstruct with the new SVs
        SVs(:,fn:end)=0;
        Xm = U*SVs*D';
        %fill missing values 
        Xfilled(missing_ind)=Xm(missing_ind);
        
        % check for convergence 
        SS = sum(sum(Xm(missing_ind).^2));
        f = abs((SS-SSold)/SSold);
        it=it+1;
    end 
    disp('f')
    disp(f)
    disp('iter')
    disp(it)
    %metrics 
     % solve for the correct number of factors 
    [~,~,D,St,X_pred,~,~] = missing_svd(Xs,fn,1,1,1e-3,1);
    Xm = St*D';
    mse = (sum((Xm(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
    wmse = find_wmse(Xs(filled_ind), Xm(filled_ind), length(filled_ind));
    Cyt = corrcoef(X_pred(filled_ind),Xm(filled_ind));
    R2=(Cyt(2,1));
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

function [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd(X,fn,center,scale,conv,max_iter)
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
    % Output 
    % S,V,D from X = SVD'
    % St = SV
    % X_pred = filled X with new values for the missing entries
    [m,n]=size(X);
    missing_ind = find(isnan(X));
    if any(isnan(X)) % there is missing data 
        Xfilled = fill_data(X); 
        SS = sum(sum(Xfilled(missing_ind).^2));
        
        f=2*conv;
        iter = 0;
        while iter<max_iter && f>conv
            iter = iter+1;
            SSold = SS;            
            % preprocess = scale and then center 
            mx = mean(Xfilled);
            if scale ==1 
                sj = sqrt(sum(sum((Xfilled).^2)));
                Xfilled = Xfilled/sj;
            end 
            if center ==1
                mx = mean(Xfilled);
                Xc = Xfilled-ones(m,1)*mx; %centering of the data done -> for PCA 
            else 
                Xc=Xfilled;
            end 
            [S,V,D]=svd(Xc);
            AllSVs = V;
            St=S*V;
             St = St(:,1:fn);
             V=V(:,1:fn);
             S=S(:,1:fn);
             D=D(:,1:fn);
            % post process = uncenter, unscale  
            if center ==1
                X_pred = St*D'+ones(m,1)*mx;
            else 
                X_pred = St*D';
            end 
            if scale ==1
                X_pred = X_pred*sj;
            end
            
            % fill misssing values 
            Xfilled(missing_ind)=X_pred(missing_ind);
            SS = sum(sum(Xfilled(missing_ind).^2));
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
        [S,V,D]=svd(Xc);
        St=S*V;
        AllSVs=V;
        St = St(:,1:fn);
        V=V(:,1:fn);
        S=S(:,1:fn);
        D=D(:,1:fn);
         
        if center ==1
            X_pred = St*D'+ones(m,1)*mx;
        else 
            X_pred = St*D';
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




