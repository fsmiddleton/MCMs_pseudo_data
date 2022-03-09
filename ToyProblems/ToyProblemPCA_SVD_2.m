%% Toy problem to create code for SVD used for matrix completion 
% Francesca Middleton, 2022-03-02

%% Import or create data matrix for testing 
clc
clear

import = 1; 
% import = 0 to create a low rank matrix, 
% import =1 to fetch an already created low rank matrix
% import =2 to import a full test matrix of spectroscopy data, 
% import =3 to import excess enthalpy data,
filename = 'RandomMatrixNoise2.xlsx';
export = 0; % either export data(1) or don't, only used for creating a low rank matrix
noise=1; %noise =1 to add noise to the creation of the matrix
rank_mat = 0; % to be assigned later 

if import ==0
    % create array
    % specify size and rank of array, choosing random mu and sigma to create
    % singular values from, and the noise  
    n=50;
    m=50;
    rank_mat = 10;
    mu = 20;
    sigma = 3;
    [Xs,Xnoise, rankU, rankV, noiseMat]=create_matrix(n,m,rank_mat, mu, sigma);
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
    if noise == 1
        T2 = readtable(filename, 'Sheet', 4);
    else
        T2 = readtable(filename, 'Sheet', 2);
    end 
    X = table2array(T2);
    T3 = readtable(filename, 'Sheet', 1);
    rank_mat = T3(1,2);
    
elseif import ==2
    % correlated matrix of spectral data, manuscript available (highly correlated data):
    % https://www.kaggle.com/sergioalejandrod/raman-spectroscopy
    T1 = readtable('raman_mix1_spectrum.xlsx', 'Sheet', 'mix_1'); 
    X = table2array(T1);%only use T1, small enough to manage

    
elseif import ==3
    % Import excess enthalpy data 
    T1 = readtable('HEmatrix298.15.2.xlsx', 'Sheet', 3); 
    X = table2array(T1);%only use T1, small enough to manage
    
else
    disp('Please choose 0, 1, 2 to use data')
end 

%%
%remove NaN row or column 
row_nan = find(all(isnan(X),2));
col_nan = find(all(isnan(X),1));
if row_nan
    X(row_nan,:)=[];
end 
if col_nan 
    X(:,col_nan)=[];
end 
n = size(X,1);
m = size(X,2);
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
remove_ind = 1:(n*m);
 
%Time to choose how to find the best rank
Gavish =0;
 % Choose sparsities used for finding the rank 
sparsities = 0.5:0.01:0.6;
%ranks to try 
fns = 1:1:20;
        
%intialise the metrics to analyse each sparsity and its final rank found 
minmse = zeros(size(sparsities,2),1);
minwmse = zeros(size(sparsities,2),1);
min_fn = zeros(size(sparsities,2),1);
X_pred_best = zeros(n,m,size(sparsities,2)); % and the X predictions
R2 = zeros(size(sparsities,2),1);
% initialise error vars 
mse = zeros(size(fns,2),size(sparsities,2));
smse = zeros(size(fns,2),size(sparsities,2));
wmse = zeros(size(fns,2),size(sparsities,2));
R = zeros(size(fns,2),size(sparsities,2));
        
j=0;
for sparsity = sparsities
    j=j+1; % for sparsity

    % Create the sparse matrix 
    if remove == 1
        [Xs,missing_ind, filled_linear_ind] = remove_matrix(X,sparsity);
    elseif remove == 2
        [Xs,missing_ind,filled_linear_ind]=fill_matrix(X,sparsity);
    elseif import ==1
        % import the relevant sparse matrix from the spreadsheet
        Ts = readtable(filename, 'Sheet', num2str(sparsity));
        Xs = table2array(Ts);
        missing_ind = find(isnan(Xs));
        filled_linear_ind = find(~isnan(Xs));
    elseif import == 3 || import ==2 
        % this case is for real experimental data 
        Xs = (X); % no removing of data necessary 
        missing_ind = find(isnan(X));
        filled_linear_ind = find(~isnan(X));
    end 
    
    % Find the optimal rank for this matrix
    
    if Gavish==1
        %Hard thresholding 
        if n/m ==1
            omega = 2.858; % omega(beta)=2.858 for n*n square matrix
        elseif n/m < 1
            omega = optimal_SVHT_coef_sigma_unknown(n/m);
        else
            beta = m/n;
            %omega = optimal_SVHT_coef_sigma_unknown(m/n); 
            omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
        end 
        [U,D,V,X_pred]=missing_svd(Xs,m,1,1e-4,200); % PCA done once, the matrix needs to be filled to do this
        
        y_med = median(diag(D)); %ymed = median singular value of the noisy matrix  
        cutoff = omega*y_med; %cutoff= tau= omega(beta)*ymed; matrix
         % Keep modes w/ sig > cutoff; rank chosen as hard cutoff
        d_original = diag(D);
        d=diag(D);
        fn = length(find(diag(D)>cutoff));
        %disp(fn)
        d(1+fn:end)=0;
        %fn = length(find(diag(D)>cutoff));
        %disp(fn)
        X_pred = U*diag(d)*V';
        SRSS = sqrt((sum((X_pred(missing_ind)-X(missing_ind)).^2)));
        minmse(j) = (sum((X_pred(missing_ind)-X(missing_ind)).^2))/length(remove_ind);
        minwmse(j)= find_wmse(X(missing_ind), X_pred(missing_ind), length(missing_ind));
        min_fn(j) = fn;
        X_pred_best(:,:,j) = X_pred;
        Cyt = corrcoef(X_pred(missing_ind),X(missing_ind));
        R2(j)=(Cyt(2,1));
        %End Gavish method 
    else
        % Iterative PCA with wMSE or MSE used to find rank 
        
        %choose which error measure to use for choosing the best PC
        winsorized_mse =1; %1=use wmse 

        
        %Find rank by minimising the mse or wmse 
        i=0;
        for fn=fns
            [U,D,V,X_pred]=missing_svd(Xs,fn,1,1e-3,1000);
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
    X_pred_boot = zeros(n,m,size(filled_linear_ind,1));
    boot_removed_col = zeros(size(filled_linear_ind,1),1);
    boot_removed_row = zeros(size(filled_linear_ind,1),1);
    k=0;% for the bootstrapping to find error bars 
    for filled_ind = filled_linear_ind' %filled_linear_ind must be a row vector for a for loop
        %Find error bars of the predictions 
        % remove a point from Xs 
        X_b = Xs;
        col = mod(filled_ind,m);
        if col ==0
            col=m;% mod(any integer*m,m)=0, but the column should be column m
        end 
        row = ceil(filled_ind/m);
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
%% Create table with results 
Results = table(["Sparsity"; (sparsities')],[ "SVs";min_fn],[ "MSE";minmse],[ "wMSE";minwmse], ["R2";R2]);
disp(Results)
%% Plots of the singular values 
clf
m=40;
sparsity = 0.5500;
if import ==1
    % import the relevant sparse matrix from the spreadsheet
    Ts = readtable(filename, 'Sheet', num2str(sparsity));
    Xs = table2array(Ts);
    missing_ind = find(isnan(Xs));
    filled_linear_ind = find(~isnan(Xs));
    disp(1)
end 
ind = find(sparsities == sparsity);

[U,D,V,X_pred_plot]=missing_svd(Xs,m,1,1e-3,1000); %uncomment to
mseplot = (sum((X_pred(missing_ind)-X(missing_ind)).^2))/length(remove_ind);
wmseplot= find_wmse(X(missing_ind), X_pred_plot(missing_ind), length(missing_ind));
smseplot = sqrt(mse);
Cyt = corrcoef(X_pred(missing_ind),Xs(missing_ind));
R=sqrt(Cyt(2,1));
%choose rank
subplot(2,1,1)
y = diag(D);
semilogy(1:m, diag(D))
xlabel('Rank')
ylabel('Singular values')

if (rank_mat)>0
    hold on 
    semilogy(rank_mat, y(rank_mat), 'ko','MarkerSize',5, 'LineWidth',5)
    legend('SVs', 'True rank')
    hold off
end 
subplot(2,1,2)
plot(1:m,cumsum(diag(D))/sum(diag(D)))
xlabel('Rank')
ylabel('Cumulative contribution to error')
if rank_mat>0
    hold on
    y2 = y(1:rank_mat);
    realranky = cumsum(y(1:m));
    plot(min_fn(ind), realranky(min_fn(ind)-fns(1)+1)/sum(diag(D)), 'r*', 'MarkerSize',15, 'LineWidth',1.5)
    legend('Cumulative error', 'Predicted rank')
    hold off
end 
sgtitle(strcat('Sparsity = ',num2str(sparsity)))
%% Plots of errors 

clf
plotno = length(sparsities);
i=0;
for sparsity = sparsities
    i=i+1;
    ind = find(sparsities == sparsity);
    subplot(4,3,i)
    titlestr = strcat('Sparsity = ',num2str(sparsity));
    
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
    plot(1:m,d_original)
    xlabel('Rank')
    ylabel('Singular value')
    xlabel('Component number')
    hold on 
    plot(rank_mat, d_original(rank_mat), 'ko')
    plot(fn, d_original(fn), 'ro')
    legend('Singular values', 'True rank', 'Predicted rank')
    
    subplot(2,1,2)
    plot(1:m,cumsum(d_original)/sum(d_original))
    xlabel('Component number')
    ylabel('Cumulative contribution to error')
end 


%% Reorder data to test the change of the order of the data (in rows and columns) and how this affects the final predictions 


%% Functions 
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

function omega = optimal_SVHT_coef_sigma_unknown(beta)
%Gavish and Donaho; finding coeff for optimal cutoff for a non-square
%matrix
% Beta = n/m = aspect ratio 
    warning('off','MATLAB:quadl:MinStepSize')
    assert(all(beta>0));
    assert(all(beta<=1));
    assert(prod(size(beta)) == length(beta)); % beta must be a vector
    
    coef = optimal_SVHT_coef_sigma_known(beta);

    MPmedian = zeros(size(beta));
    for i=1:length(beta)
        MPmedian(i) = MedianMarcenkoPastur(beta(i));
    end

    omega = coef ./ sqrt(MPmedian);
end

function lambda_star = optimal_SVHT_coef_sigma_known(beta)
    assert(all(beta>0));
    assert(all(beta<=1));
    assert(prod(size(beta)) == length(beta)); % beta must be a vector
    
    w = (8 * beta) ./ (beta + 1 + sqrt(beta.^2 + 14 * beta +1)); 
    lambda_star = sqrt(2 * (beta + 1) + w);
end
function med = MedianMarcenkoPastur(beta)
    MarPas = @(x) 1-incMarPas(x,beta,0);
    lobnd = (1 - sqrt(beta))^2;
    hibnd = (1 + sqrt(beta))^2;
    change = 1;
    while change & (hibnd - lobnd > .001),
      change = 0;
      x = linspace(lobnd,hibnd,5);
      for i=1:length(x),
          y(i) = MarPas(x(i));
      end
      if any(y < 0.5),
         lobnd = max(x(y < 0.5));
         change = 1;
      end
      if any(y > 0.5),
         hibnd = min(x(y > 0.5));
         change = 1;
      end
    end
    med = (hibnd+lobnd)./2;
end

function I = incMarPas(x0,beta,gamma)
    if beta > 1,
        error('betaBeyond');
    end
    topSpec = (1 + sqrt(beta))^2;
    botSpec = (1 - sqrt(beta))^2;
    MarPas = @(x) IfElse((topSpec-x).*(x-botSpec) >0, ...
                         sqrt((topSpec-x).*(x-botSpec))./(beta.* x)./(2 .* pi), ...
                         0);
    if gamma ~= 0,
       fun = @(x) (x.^gamma .* MarPas(x));
    else
       fun = @(x) MarPas(x);
    end
    I = quadl(fun,x0,topSpec);
    
    function y=IfElse(Q,point,counterPoint)
        y = point;
        if any(~Q),
            if length(counterPoint) == 1,
                counterPoint = ones(size(Q)).*counterPoint;
            end
            y(~Q) = counterPoint(~Q);
        end
        
    end
end


