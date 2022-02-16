% correlated matrix of spectral data, manuscript available (highly correlated data):
% https://www.kaggle.com/sergioalejandrod/raman-spectroscopy


%% Import or create data matrix for testing 
clc
clear

import = 0;
if import ==1
    T1 = readtable('raman_mix1_spectrum.xlsx', 'Sheet', 'mix_1');
    T2 = readtable('raman_mix2_spectrum.xlsx', 'Sheet', 'mix_2');
    T3 = readtable('raman_mix3_spectrum.xlsx', 'Sheet', 'mix_3');
    T4 = readtable('raman_mix4_spectrum.xlsx', 'Sheet', 'mix_4');

    T = cat(1, T1, T2, T3, T4); % concatenate arrays vertically 
    X = table2array(T1);
    n = size(X,2);
    m = size(X,1);
else 
    %create array
    %speficy size and rank of array, choosing random mu and sigma to create
    %singular values from 
    n=200;
    m=150;
    rank = 7;
    mu = 10;
    sigma = 5;
    [X, rankU, rankV]=create_matrix(n,m,rank, mu, sigma);
    
end 

% remove data or fill a matrix 
remove = 1;
remove_ind = 1:(n*m);

%choose how sparse test matrices will be and the number of PCs tested 
sparsities = [0.1,0.2,0.3,0.4,0.5,0.6];
fns = [1,2,3,4,5,6,7,8,9,10];

%choose which error measure to use for choosing the best PC
winsorized_mse =1; %1=use wmse 

% initialise error vars 
mse = zeros(size(fns,2),1);
smse = zeros(size(fns,2),1);
wmse = zeros(size(fns,2),1);
minmse = zeros(size(sparsities,2),1);
minwmse = zeros(size(sparsities,2),1);
min_fn = zeros(size(sparsities,2),1);
X_pred_best = zeros(n,m,size(sparsities,2)); % and the X predictions

j=0;
for sparsity=sparsities
    j=j+1;
    if remove ==1
        [Xs,missing_ind] = remove_matrix(X,sparsity);
    else 
        [Xs,missing_ind]=fill_matrix(X,sparsity);
    end 
    
    i=0;
    for fn=fns
        [U,D,V,X_pred]=missing_svd(Xs,fn,1,1e-3,1000);
        i=i+1;
        SRSS = sqrt((sum((X_pred(missing_ind)-X(missing_ind)).^2)));
        mse(i) = (sum((X_pred(missing_ind)-X(missing_ind)).^2))/length(remove_ind);
        wmse(i)= find_wmse(X(missing_ind), X_pred(missing_ind), length(missing_ind));
        smse(i) = sqrt(mse(i));
    end
    
    minmse(j) = min(mse);
    minwmse(j)=min(wmse);
    if winsorized_mse ==1
        min_index = find(wmse==minwmse(j));
    else
        min_index = find(mse==minmse(j));
    end 
    min_fn(j) = fns(min_index);
    [U,D,V,X_pred_best(:,:,j)]=missing_svd(Xs,min_fn,1,1e-3,1000);
end 

 


%% Plots of the singular values 
subplot(2,1,1)
semilogy(1:n, diag(D))
subplot(2,1,2)
plot(1:n,cumsum(diag(D))/sum(diag(D)))
%% Functions 
function [X, rankU, rankV]=create_matrix(n,m,r, mu, sigma)
    % create a matrix of size nxm with rank =r using SVD 
    
    sing_vals = normrnd(mu, sigma, r,1); % generate r singular values that are normally distributed 
    sing_vals = sort(sing_vals, 'descend'); % sorted from highest to lowest
    D = zeros(n,m);%initialise D 
    for i = 1:r
        disp(i)
        D(i,i) = sing_vals(i); %populate array with values 
    end 
    % fill the matrices U,V with random entries
    U = rand(n,n);
    V = rand(m,m);
    % find orthonormal bases of U and V 
    rankU = rank(U);
    rankV = rank(V);
    U = orth(U);
    V = orth(V);
    %fill the final matrix 
    X = U*D*V'; %create the matrix of rank r
end 

function [X_sparse,remove_ind]=fill_matrix(X, perc_remove)
    % fill data into a matrix until you have reached 1-perc_remove, save all
    % indices not filled 
    [n,m]=size(X);
    num_removed = floor(m*n*perc_remove);
    %create NaN matrix 
    X_sparse = NaN(n,m);
    X_filled_indices = zeros(n,m);
    % fill in one entry in every row and column, randomly 
    filled_cols = [];
    cols=randperm(m); % sort the columns in a random order 
    for i = 1:n
        for j =1:cols 
            X_sparse(i, j)=X(i,j);
            X_filled_indices(i,j)=1;
        end 
    end 
    
    for k= 1:(m*n-num_removed) % remove this many values from array X 
        i = randi(m,1); %matrix inidices start at 1 in MATLAB
        j = randi(n,1);
        while X_filled_indices(i,j)==1
            %repeat the filling if this index is already filled 
            i = randi(n,1); %row
            j = randi(m,1);%col
        end 
        X_sparse(i,j)= X(i,j); %reassign value as original X value 
    end 
    remove_ind = find(isnan(X_sparse));
end 

function [X_sparse,remove_ind]=remove_matrix(X,perc_remove)
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
end 

function [X_filled]=fill_data(X)
    [m,n]=size(X);
    missing_ind = (isnan(X));
    
    X_filled=X;
    X_filled(missing_ind)=0; %fill NaN values with 0 
    
    [i, j]=find(isnan(X));% returns rows and columns with nonzero elements 
    
    mean_col = sum(X_filled,1)./(ones(1,n)*m-sum(missing_ind,1)); %columns are dimension 1
    mean_row= sum(X_filled,2)./(ones(1,m)*n-sum(missing_ind,2)); % rows are dimension 2 
    for k =1:length(i) % for all NaN elements that exist, loop through them to replace with means 
        X_filled(i(k),j(k))=(mean_row(i(k))+mean_col(j(k)))/2;
    end 
end 

function [S,V,D,X_pred]=missing_svd(X,fn,center,conv,max_iter)
    [m,n]=size(X);
    miss_ind = find(isnan(X));
    if any(isnan(X)) % there is missing data 
        Xf = fill_data(X); 
        mx = mean(Xf);
        SS = sum(sum(Xf(miss_ind).^2));
        
        f=2*conv;
        iter = 1;
        while iter<max_iter && f>conv
            SSold = SS;
            if center ==1
                mx = mean(Xf);
                %Xc = normalize(Xf); %normalizing 
                Xc = Xf-ones(m,1)*mx; %centering
            else 
                Xc=Xf;
            end 
            [S,V,D]=svd(Xc);
            S=S*V;
            S=S(:,1:fn);
            D=D(:,1:fn);
            if center ==1
                X_pred = S*D'+ones(m,1)*mx;
            else 
                X_pred = S*D';
            end 
            Xf(miss_ind)=X_pred(miss_ind);
            SS = sum(sum(Xf(miss_ind).^2));
            f = abs(SS-SSold)/(SSold);
            iter = iter+1;
        end 
    else 
        Xf=X;
        if center ==1
            mx = mean(Xf);
            %Xc = normalize(Xf); %normalizing 
            Xc = Xf-ones(m,1)*mx; %centering
        else 
            Xc=X;
        end 
        [S,V,D]=svd(Xc);
        S=S*V;
        S=S(:,1:fn);
        D=D(:,1:fn);
        if center ==1
            X_pred = S*D'+ones(m,1)*mx;
        else 
            X_pred = S*D';
        end 
    end %end if  else 
end % end function 

function [wmse]=find_wmse(X,X_pred,no_missing)
% find the winsorized mse for the matrices 
    % turn 
    % find outliers
    perc5 = prctile(X_pred,5, 'all');
    perc95 = prctile(X_pred, 95, 'all');
    %reassign
    X_pred(X_pred<perc5)=perc5;
    X_pred(X_pred> perc95)=perc95;
    wmse = (sum((X_pred-X).^2))/no_missing;
end 