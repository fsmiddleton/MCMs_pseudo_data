%% PCA done using SVD 
% Francesca Middleton 2022-04-25 

%% Make data 

% Parameters for the 3-way array 
dim1 = 4;  % Size of the first dimension
dim2 = 5;  % Size of the second dimension
dim3 = 6;  % Size of the third dimension
dim = [dim1,dim2,dim3];
missing_percentage = 0.2;  % Percentage of missing entries

% Generate a random 3-way array
X = rand(dim1, dim2, dim3);

% Determine the number of missing entries
num_elements = numel(X(:,:,1));
num_missing = round(missing_percentage * num_elements);

% Randomly select indices to set as missing (NaN)
for i =1:dim3
    missing_indices(:,i) = randperm(num_elements, num_missing);
    Xtemp = X(:,:,i);
    Xtemp(missing_indices(:,i))=NaN;
    X_missing(:,:,i) = Xtemp;
end
% Set the selected indices to NaN


% Display the result
disp('Random 3-way array with missing entries:');
disp(X_missing);

 
%% Matrix completion step 
% Set some variables 
fns = [3:5]; % ranks to test: must be < dim2
n = 10; % number of iterations 
scale = 1; % choose whether to scale
center = 1; % choose whether to center 

% Intialise the metrics to analyse each sparsity and its final rank found 
% vectors, one for each matrix that is sliced 
minmse = zeros(size(dim3,2),1); 
minwmse = zeros(size(dim3,2),1);
min_fn = zeros(size(dim3,2),1);
X_pred_best = zeros(dim1,dim2,size(dim3,2)); % and the X predictions
R2 = zeros(size(dim3,2),1);
% initialise error vars 
% column for each matrix (fn vs sparsity)
mse = zeros(size(fns,2),size(dim3,2));
smse = zeros(size(fns,2),size(dim3,2));
wmse = zeros(size(fns,2),size(dim3,2));
R = zeros(size(fns,2),size(dim3,2));
tiledlayout(dim3,1)
i=0; % counter 

for dimension3 = 1:dim3 
    i = i+1;
    X_truth = X(:,:,i);
    X_temp = X_missing(:,:,i);
    missing_ind = missing_indices(:,i);
    Xs=X_truth(missing_ind);
    % complete matrix
    % Iterative PCA with wMSE or MSE used to find rank 

    %choose which error measure to use for choosing the best PC
    winsorized_mse =1; %1=use wmse 

    %Find rank by minimising the mse or wmse 
    j=0;
    for fn=fns
        j=j+1;
        [U,D,V,X_pred]=missing_svd(X_temp,fn,1e-3,n,scale,center);
        X_pred_all(:,:,i,j)= X_pred;
        
        SRSS =sum(sum((X_pred(missing_ind)-Xs).^2))/length(missing_ind);
        % disp('next dim and rank')
        % disp(i)
        % disp(j)
        % disp(SRSS)
        mse(i,j) = sum(sum((X_pred(missing_ind)-Xs).^2))/length(missing_ind);
        wmse(i,j)= find_wmse(Xs, X_pred(missing_ind), length(missing_ind));
        smse(i,j) = sqrt(mse(i));
    end
    minmse(i) = min(mse(i,:));
    minwmse(i)=min(wmse(i,:));

    min_index = find(wmse(i,:)==minwmse(i));
    if length(min_index)>1
       min_index=1;
    end 
    min_fn(i) = fns(min_index);
    [U,D,V,X_pred_best(:,:,i)]=missing_svd(X_temp,min_fn(i),1e-3,n,scale,center);
    y = diag(D);
    nexttile
    semilogy(fns, diag(D))
    xlabel('Rank')
    ylabel('Singular values')
    xticks(fns)
end 

%% Functions 

function [X_filled]=fill_data(X)
    [m,n]=size(X);
    missing_ind = (isnan(X));
    [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
    X_filled=X;
    X_filled(missing_ind)=0; %fill NaN values with 0
    mean_col = sum(X_filled,1)./(sum(missing_ind,1)); %columns are dimension 1
    mean_row= sum(X_filled,2)./(sum(missing_ind,2)); % rows are dimension 2 
    for k =1:length(i) % for all NaN elements that exist, loop through them to replace with means 
        X_filled(i(k),j(k))=(mean_row(i(k))+mean_col(j(k)))/2;
    end 
end 

function [S,V,D,X_pred]=missing_svd(X,fn,conv,max_iter, scale, center)
    [m,n]=size(X);
    miss_ind = find(isnan(X));
    if length(miss_ind)>0 % there is missing data 
        Xf = fill_data(X); 
        mx = mean(Xf);
        SS = sum(sum(Xf(miss_ind).^2));
        
        f=2*conv;
        iter = 1;
        while iter<max_iter && f>conv
            SSold = SS;
            if center ==1
                mx = mean(Xf);
                Xc = Xf-ones(m,1)*mx; %centering of the data done -> for PCA 
            else 
                Xc=Xf;
            end 
            [S,V,D]=svd(Xc);
            S=S*V;
            V=V(:,1:fn);
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


