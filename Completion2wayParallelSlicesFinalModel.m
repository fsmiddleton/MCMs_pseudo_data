%% 2-way arrays 
% Francesca Middleton, 2022-03-02


%% Import .m file with 3-way array
clc
clear
% decide on which interval and temperature to evaluate 
interval = 0.05;
T =298.15;
% choose which array you want to work with 
filename = strcat('HEData3wayPolySmall-',num2str(interval), '-', num2str(T), '.mat');
%filename = strcat('HEData3wayPolyAll-',num2str(interval), '-', num2str(T), '.mat');
load(filename)
conc_interval = interval:interval:(1-interval);
indend = 20;
X = HE_data_sparse(1:indend,1:indend,1:length(conc_interval));
Xsign = sign(X);
Xscale = Xsign.*log(Xsign.*X);
dim = size(X);
dim1 = dim(1);
dim2 = dim(2);
dim3 = dim(3);
percmiss = length(find(isnan(X)))/(dim1*dim2*dim3)*100;
percobs = length(find(~isnan(X)))/(dim1*dim2*dim3)*100;

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
%% fill the X array with previously created answers 
load('2waySVD-50comps-scale1-threshold60par-moreiter-huberloss1-LOOCV-Xscale-T=298.15-fillmethod=rec-02-Nov-2022.mat', 'X_pred')
%fill Xscale values with those created 
indend = size(X,1);
%% Rank determined from SVs
% Import other data for most missing to allow variables to be declared
concentrations = conc_interval;
%filename = 'TestSMALLHEMatrix13June298.15.xlsx';
fn = min([dim1,dim2]);
y1 = zeros(fn,length(concentrations));
y2 = zeros(fn,length(concentrations));

%necessary for plots
indexplot=20;
whichX = 'sign';
if strcmp(whichX,'scale')
    Xss=Xscale;
else 
    Xss=Xsign;
end
i=0;
for c = (concentrations)
    i=i+1;
    %Import data 
    %T = readtable(filename, 'Sheet', num2str(c));
    Xs=Xss(:,:,i);
    %SVD
    Xf = filldata3(Xs, 'uni',mixtures,c,whichX,T);
    [~,D,~,~,~, ~]=missing_svd(Xf,fn,1,1,1e-8,1,0);
    %Plot
    y1(:,i) = diag(D);
    y2(:,i) = cumsum(y1(:,i))/sum(y1(:,i));
    subplot(2,1,1)
    plot(1:indexplot, y1(1:indexplot,i))
    subplot(2,1,2)
    plot(1:indexplot,y2(1:indexplot,i))
    hold on
end 

hold off 
filenamescree = strcat('Scree-All-',whichX,'-',num2str(T),'K-interval=',num2str(interval),'.mat');
save(filenamescree);
%% LOOCV to find rank using mse  
%time the code 
tic
% declare ranks to be tested 
fns =[8:16];
concentrations=conc_interval;
Xscale = log(sign(X).*(X)).*sign(X);
Xsign = sign(X);
Xs = Xscale;
tempX = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
[row,col] = find(~isnan(tempX));
filled_ind = find(~isnan(tempX));

% vars for missing_parafac3
maxiter = 20000;
scale = 1;
center = 1;
conv = 1e-10;
winsorized_mse = 0;
fillmethod = 'reg';
orth = 1;
whichX = 'scale';

% declare vars for analysis 
mse_LOOCV = zeros(length(concentrations),length(fns));
wmse_LOOCV = zeros(length(concentrations),length(fns));
RAD_LOOCV = zeros(length(concentrations),length(fns)); % relative absolute deviation 
min_mse_LOOCV = zeros(length(concentrations),1);
min_wmse_LOOCV = zeros(length(concentrations),1);
X_pred_LOOCV = zeros(dim1,dim2,length(filled_ind));
LOOCV_removed_col = zeros(length(filled_ind),1);
LOOCV_removed_row = zeros(length(filled_ind),1);
Xm_boot=zeros(length(concentrations), length(fns), length(filled_ind));
Xm_boot2=zeros(length(concentrations), length(fns), length(filled_ind));

conc = concentrations;
 if strcmp(whichX, 'scale') 
    Xs = Xscale;
 elseif strcmp(whichX,'sign') 
    Xs = Xsign;
 else
    Xs = X;
end 

%declare vars with size dependent on the array used 
%loop through ranks
disp('fn')
fnind=0;
for fn = fns 
    disp(fn)
    fnind = fnind + 1; 
        
    X_b = Xs;

    %perform iterative PCA on the slightly more empty matrix 
    [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd_par(X_b,fn,center,scale,conv,maxiter, 1,fillmethod,mixtures,whichX,conc,T);
    filenamesave = strcat('2waySVD-final-50comps-threshold60par-X',whichX,'-T=',num2str(T), '-fillmethod=',fillmethod,'-fn=',num2str(fn),'-',  date, '.mat');
    save(filenamesave)

end % END FN
    




%% Functions 

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
    errors = errors(:);
    perc5=prctile(errors,5,'all');
    perc95=prctile(errors,95,'all');
    errors(errors<perc5)=perc5;
    errors(errors>perc95)=perc95;
    wmse = (sum((errors).^2))/count;
end 

function [i,j,k]=findnan3(X)
% Input 
% X = data array 
% Output 
% i, j, k = indexes of nan values in X

% isnan(X) returns logical array 
% all(isnan(X),1) is also a logical array (1x30x40) - each row
% find(all(isnan(X),1)) returns linear indices -> reshape array to correct
% dimension to find missing slabs 
% findnan3 finds the nan indices of a 3-way array 
    dim = size(X);
    dim1 = dim(1);
    dim2 = dim(2);
    dim3 = dim(3);

    i=[];
    j=[];
    k=[];
    for d = 1:dim3
        % done per z slice 
        Xtemp = X(:,:,d);
        [itemp, jtemp]= find(isnan(reshape(Xtemp,[dim1,dim2])));
        i = [i; itemp];
        j = [j;jtemp];
        ktemp = ones(length(itemp),1)*d;
        k = [k;ktemp];
    end 
end

