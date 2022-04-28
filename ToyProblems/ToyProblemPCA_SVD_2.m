%% Toy problem to create code for SVD used for matrix completion 
% Francesca Middleton, 2022-03-02

%% Import or create data matrix for testing 
clc
clear

import = 1; 
% import = 0 to create a low rank matrix, 
% import =1 to fetch an already created low rank matrix (2-way) 
% import =2 to import excess enthalpy data (2-way),
% import = 3 to import 3-way data 
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
            [Xsparse,missing_ind,filled_linear_ind]=fill_matrix(X,i);
            Tsparse = array2table(Xsparse);
            writetable(Tsparse, filename, 'Sheet', num2str(i))
        end 
    end
    
        
elseif import ==1
    %import the data from the previously created random matrix 
    filename = 'RandomMatrixNoNoise.xlsx'; % choose file here
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
    
elseif import ==3
    % Import multiple sheets, 3 way array 
    filename = 'ToyProblemData3DFull.xlsx'; % created using parafac code 
    % find sheet names 
    sheets = sheetnames(filename);
    % iterate for the number of sheets
    dim3 =0;
    for s=1:length(sheets)
        [sNum,tf]=str2num(s);
        % if the name of the sheet was a number, it can be converted and
        % contains data for the problem 
        if tf == 1
            T = readtable(filename, 'Sheet', sNum);
            X(:, :, sNum)=table2array(T); 
        end 
    end 
    dim = size(X);
    
else
    disp('Please choose 0, 1, 2 to use data')
end 

%% remove NaN rows or columns
% saves the columns and rows with Nan values 
row_nan = find(all(isnan(X),2));
col_nan = find(all(isnan(X),1));
if row_nan
    X(row_nan,:)=[];
end 
if col_nan 
    X(:,col_nan)=[];
end 
dim1 = dim(1);
dim2 = dim(2);
rank_mat=table2array(rank_mat);

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
%%
% remove data or fill a matrix 
%filename specified in first block of code for import ==1 
remove = 0; % 1=remove to create sparse matrix, 2 = fill a sparse matrix, 0 = no change to matrix
remove_ind = 1:(dim1*dim2);
 
%Time to choose how to find the best rank
Gavish =0;
 % Choose sparsities used for finding the rank 
intervals = 0.5:0.01:0.6;
%ranks to try 
fns = 1:1:20;
%choose whether to reorder matrix or not, 1 = true 
reorder = 0;
        
%intialise the metrics to analyse each sparsity and its final rank found 
% vectors, one for each matrix that is sliced 
minmse = zeros(size(intervals,2),1); 
minwmse = zeros(size(intervals,2),1);
min_fn = zeros(size(intervals,2),1);
X_pred_best = zeros(dim1,dim2,size(intervals,2)); % and the X predictions
R2 = zeros(size(intervals,2),1);
% initialise error vars 
% column for each matrix (fn vs sparsity)
mse = zeros(size(fns,2),size(intervals,2));
smse = zeros(size(fns,2),size(intervals,2));
wmse = zeros(size(fns,2),size(intervals,2));
R = zeros(size(fns,2),size(intervals,2));



j=0; % counter for intervals

for composition = intervals
    j=j+1; % for sparsity

    % Create the sparse matrix for toy problem 
    if remove == 1
        [Xs,missing_ind, filled_linear_ind] = remove_matrix(X,composition);
    elseif remove == 2
        [Xs,missing_ind,filled_linear_ind]=fill_matrix(X,composition);
    elseif remove ==0 && import ==1
        % import the relevant sparse matrix from the spreadsheet
        Ts = readtable(filename, 'Sheet', num2str(composition));
        Xs = table2array(Ts);
        if reorder ==1
            % redorder matrix randomly, choose how to here 
            p = reorder_He(Xs, "cols");
            Xs(:) = Xs(p);
            X(:)=X(p);% reorder original matrix to allow the mse, wmse and R2 to work 
        end 
        missing_ind = find(isnan(Xs));
        filled_linear_ind = find(~isnan(Xs));
    elseif import == 3 || import ==2 || remove ==0
        % this case is for real experimental data 
        Xs = (X); % no removing or filling of data necessary 
        missing_ind = find(isnan(X));
        filled_linear_ind = find(~isnan(X));
    end
    
    
    % Find the optimal rank for this matrix
    
    if Gavish==1
        [min_fn(j),minmse(j),minwmse(j),R2(j), X_pred_best(:,:,j)] = solveGavish(Xs, dim1, dim2);
         
    else
        % Iterative PCA with wMSE or MSE used to find rank 
        
        %choose which error measure to use for choosing the best PC
        winsorized_mse =1; %1=use wmse 

        %Find rank by minimising the mse or wmse 
        i=0;
        for fn=fns
            %[U,D,V,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
            [U,D,V,X_pred]=missing_svd(Xs,fn,1,1,1e-3,1000);
            i=i+1;
            SRSS = sqrt((sum((X_pred(missing_ind)-X(missing_ind)).^2)));
            mse(i,j) = (sum((X_pred(missing_ind)-X(missing_ind)).^2))/length(remove_ind);
            wmse(i,j)= find_wmse(X(missing_ind), X_pred(missing_ind), length(missing_ind));
            smse(i,j) = sqrt(mse(i));
            Cyt = corrcoef(X_pred(missing_ind),X(missing_ind));
            R(i,j)=sqrt(Cyt(2,1));

        end
        minmse(j) = min(mse(:,j));
        minwmse(j)=min(wmse(:,j));
        
        if winsorized_mse ==1
            min_index = find(wmse(:,j)==minwmse(j));
        else
            min_index = find(mse(:,j)==minmse(j));
        end 
        min_fn(j) = fns(min_index);
        R2(j) = R(min_index,j)^2;
        [U,D,V,X_pred_best(:,:,j)]=missing_svd(Xs,min_fn(j),1,1e-3,1000);
    end
    
    %Error bars
    X_pred_boot = zeros(dim1,dim2,size(filled_linear_ind,1));
    boot_removed_col = zeros(size(filled_linear_ind,1),1);
    boot_removed_row = zeros(size(filled_linear_ind,1),1);
    k=0;% for the bootstrapping to find error bars 
    for filled_ind = filled_linear_ind' %filled_linear_ind must be a row vector for a for loop
        %Find error bars of the predictions 
        % remove a point from Xs 
        X_b = Xs;
        col = mod(filled_ind,dim2);
        if col ==0
            col=dim2;% mod(any integer*m,m)=0, but the column should be column m
        end 
        row = ceil(filled_ind/dim2);
        X_b(filled_ind) = nan;
        if find(X_b(:,col)) & find(X_b(row,:)) %ensure at least one value in each column and row
            k=k+1;
            boot_removed_col(k) = col;
            boot_removed_row(k) = row;
            %if there are other entries in the rows and columns,
            %perform iterative PCA on the slightly more empty matrix 
            X_b=fill_data(X_b);
            [U,D,V,X_pred_boot(:,:,k)]=missing_svd(X_b,min_fn,1,1e-3,1000);%iterative PCA using the known rank 
        end 
    end 
    % remove the predictions that were not made from the averaging and
    % finding variance 
    
    for l=1:k
        X_pred_best_boot(boot_removed_row(l), boot_removed_col(l)) = mean(X_pred_boot(boot_removed_row(l), boot_removed_col(l), :));
        X_var_best_boot(boot_removed_row(l), boot_removed_col(l)) = var(X_pred_boot(boot_removed_row(l), boot_removed_col(l), :));
    end 
    
end 
% Create table with results 
Results = table(["Sparsity"; (intervals')],[ "SVs";min_fn],[ "MSE";minmse],[ "wMSE";minwmse], ["R2";R2]);
disp(Results)
%% Plots of the singular values 
clf
dim2=40;
composition = 0.5200;
if import ==1
    % import the relevant sparse matrix from the spreadsheet
    Ts = readtable(filename, 'Sheet', num2str(composition));
    Xs = table2array(Ts);
    missing_ind = find(isnan(Xs));
    filled_linear_ind = find(~isnan(Xs));
    disp(1)
end 
ind = find(intervals == composition);

[U,D,V,X_pred_plot]=missing_svd(Xs,dim2,1,1e-3,1000); %uncomment to
mseplot = (sum((X_pred(missing_ind)-X(missing_ind)).^2))/length(remove_ind);
wmseplot= find_wmse(X(missing_ind), X_pred_plot(missing_ind), length(missing_ind));
smseplot = sqrt(mse);
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
sgtitle(strcat('Sparsity = ',num2str(composition)))
%% Plots of errors 
clf
plotno = length(intervals);
a = ceil(sqrt(plotno));
b = ceil(plotno/a);
i=0;
for composition = intervals
    i=i+1;
    ind = find(intervals == composition);
    subplot(a,b,i)
    titlestr = strcat('Sparsity = ',num2str(composition));
    plot(fns, wmse(:,ind))
    hold on
    plot(fns,mse(:,ind), 'm')
    plot(rank_mat,mse(rank_mat-fns(1)+1,ind), 'ko', 'MarkerSize',5, 'LineWidth',5)
    plot(min_fn(ind),wmse(min_fn(ind)-fns(1)+1,ind), 'r*', 'MarkerSize',15, 'LineWidth',1.5)
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

function  [fn,minmse,minwmse,R2, X_pred] = solveGavish(Xs, dim1, dim2)
% Gavish Hard thresholding 
        if dim1/dim2 ==1
            omega = 2.858; % omega(beta)=2.858 for n*n square matrix
        elseif dim1/dim2 < 1
            beta = dim1/dim2;
            omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
        else
            beta = dim2/dim1; 
            omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
        end 
        [U,D,V,X_pred]=missing_svd(Xs,dim2,1,1e-4,200); % PCA done once, the matrix needs to be filled to do this
        
        y_med = median(diag(D)); %ymed = median singular value of the noisy matrix  
        cutoff = omega*y_med; %cutoff= tau= omega(beta)*ymed; matrix
         % Keep modes w/ sig > cutoff; rank chosen as hard cutoff
        d_original = diag(D);
        d=diag(D);
        fn = length(find(diag(D)>cutoff));
        d(1+fn:end)=0;
        X_pred = U*diag(d)*V';
        SRSS = sqrt((sum((X_pred(missing_ind)-X(missing_ind)).^2)));
        minmse= (sum((X_pred(missing_ind)-X(missing_ind)).^2))/length(remove_ind);
        minwmse= find_wmse(X(missing_ind), X_pred(missing_ind), length(missing_ind));
        Cyt = corrcoef(X_pred(missing_ind),X(missing_ind));
        R2=(Cyt(2,1));
        %End Gavish method
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

function [X_sparse,remove_ind, fill_ind]=fill_matrix(X, perc_remove)
    % fill data into a matrix until you have reached 1-perc_remove, save all
    % indices not filled 
    [n,m]=size(X);
    num_removed = floor(m*n*perc_remove);
    %create NaN matrix 
    X_sparse = NaN(n,m);
    X_filled_indices = zeros(n,m);
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
        X_filled_indices(i,j)=1;
    end 
  
    
    for k= 1:(m*n-num_removed-m) % remove this many values from array X, 
                                  % having filled at least one entry into every row and column 
        i = randi(n,1); %matrix inidices start at 1 in MATLAB
        j = randi(m,1);
        while X_filled_indices(i,j)==1
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

function [X_filled]=fill_data(X)
    [m,n]=size(X);
    missing_ind = (isnan(X));
    [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
    X_filled=X;
    X_filled(missing_ind)=0; %fill NaN values with 0
    %disp(X_filled)
    %disp(length(i))
    mean_col = sum(X_filled,1)./(ones(1,n)*m-sum(missing_ind,1)); %columns are dimension 1
    mean_row= sum(X_filled,2)./(ones(1,m)*n-sum(missing_ind,2)); % rows are dimension 2 
    
    for k =1:length(i) % for all NaN elements that exist, loop through them to replace with means 
        X_filled(i(k),j(k))=(mean_row(i(k))+mean_col(j(k)))/2;
         
    end 
end 

function [S,V,D,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
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
            % preprocess = scale and then center 
            if scale ==1 
                sj = sqrt(sum(sum((Xf).^2)));
                Xf = Xf/sj;
            end 
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
            % post process = uncenter, unscale  
            if center ==1
                X_pred = S*D'+ones(m,1)*mx;
            else 
                X_pred = S*D';
            end 
            if scale ==1
                X_pred = X_pred*sj;
            end
            
            % fill misssing values 
            Xf(miss_ind)=X_pred(miss_ind);
            SS = sum(sum(Xf(miss_ind).^2));
            f = abs(SS-SSold)/(SSold);
            iter = iter+1;
        end 
    else 
        % no missing data 
        Xf=X;
        if scale ==1 
            sj = sqrt(sum(sum((Xf).^2)));
            Xf = Xf/sj;
        end
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
        if scale ==1
                X_pred = X_pred*sj;
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




