%% 2-way arrays 
% Francesca Middleton, 2022-03-02


%% Import .m file with 3-way array
clc
clear
% decide on which interval and temperature to evaluate 
interval = 0.05;
T =298.15;
% choose which array you want to work with 
%filename = 'HEData3wayPolyMid-0.05-298.15.mat';
filename = 'HEArrayPolySmall298.15.mat';
%filename = strcat('HEData3wayPolySmall-',num2str(interval), '-', num2str(T), '.mat');
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
load('2waySVD-final-50comps-threshold60par-Xscale-T=298.15-fillmethod=reg-fn=14-03-Nov-2022.mat', 'X_pred')
%fill Xscale values with those created 
indend = size(X_pred,1);
X(1:indend,1:indend,:) = X_pred;
%% Rank determined from SVs
% Import other data for most missing to allow variables to be declared
concentrations = conc_interval;
%filename = 'TestSMALLHEMatrix13June298.15.xlsx';
fn = min([dim1,dim2]);
y1 = zeros(fn,length(concentrations));
y2 = zeros(fn,length(concentrations));

%necessary for plots
indexplot=size(X,1);
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
%filenamescree = strcat('Scree-All-',whichX,'-',num2str(T),'K-interval=',num2str(interval),'.mat');
%save(filenamescree);
%% LOOCV to find rank using mse  
%time the code 
tic
% declare ranks to be tested 
fns =[8:2:26];
concentrations=conc_interval;
Xscale = log(sign(X).*(X)).*sign(X);
Xsign = sign(X);
Xs = Xscale;
tempX = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
[row,col] = find(~isnan(tempX));
filled_ind = find(~isnan(tempX));

% vars for missing_parafac3
maxiter = 15000;
for fillmethod = ["rec","reg"] 
scale = 1;
center = 1;
conv = 1e-10;
winsorized_mse = 0;
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
    error_LOOCV = zeros(length(filled_ind),size(X,3));
    RAD = zeros(length(filled_ind),size(X,3));
    error_LOOCV2 = zeros(length(filled_ind),size(X,3));
    RAD2 = zeros(length(filled_ind),size(X,3));
    parfor k =  1:length(filled_ind)%52 % must be a row vector for a for loop    %56,60,102
        filled_index = filled_ind(k);
        % remove a point from Xs
        X_b = Xs;
        X_b(row(k),col(k),:) = nan;
        X_b(col(k),row(k),:) = nan;
        if find(~isnan(X_b(:,col(k)))) & find(~isnan(X_b(row(k),:))) & Xs(filled_index)~=0 %ensure at least one value in each column and row

            %perform iterative PCA on the slightly more empty matrix 
            [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd_par(X_b,fn,center,scale,conv,maxiter, 1,fillmethod,mixtures,whichX,conc,T);
            %[X_pred,iters,F,err] = missing_svd(X_b,fn,maxiter,conv,scale,center,fillmethod,orth, mixtures,conc, whichX,T);
            %Factors{fn,c,k} = F;
            error_LOOCV(k,:) = Xs(row(k),col(k),:)-X_pred(row(k),col(k),:);
            error_LOOCV2(k,:) = Xs(col(k),row(k),:)-X_pred(col(k),row(k),:);
            Xm_boot(:, fnind, k) = X_pred(row(k),col(k),:);
            Xm_boot2(:, fnind, k) = X_pred(col(k),row(k),:);
            if Xs(row(k),col(k))~=0
                RAD(k,:) = error_LOOCV(k,:)./reshape(Xs(row(k),col(k),:),1,size(X,3));
                RAD2(k,:) = error_LOOCV2(k,:)./reshape(Xs(col(k),row(k),:),1,size(X,3));
            end

        end
    end
    % mse for this composition and rank 
    mse_LOOCV(fnind)= sum(sum(error_LOOCV.^2))/prod(size((find(error_LOOCV))));
    wmse_LOOCV(fnind) = find_wmse_error([error_LOOCV; error_LOOCV2], prod(size((find(error_LOOCV))))*2);
    %absolute average deviation
    RAD_LOOCV(fnind) = sum(sum(abs(RAD)))/prod(size(find(RAD)));

end % END FN
    
filenamesave = strcat('2waySVD-27comps-adding7comps-threshold60par-LOOCV-X',whichX,'-T=',num2str(T), '-fillmethod=',fillmethod,'-',  date, '.mat');
save(filenamesave)
end 
%out = Xm_boot(:,19,52);

%% Plot errors 
clf
conc=0.1;
c=conc*10;
plot(fns, (RAD_LOOCV(c,:)))

%% Find fits of the models to the data 

% declare ranks to be tested 
fns =[4:1:10, 12:2:18];
concentrations=conc_interval;
Xscale = log(sign(X).*(X)).*sign(X);
Xsign = sign(X);
%Xscale must have zeros on the diagonal too
for i =1:dim3
    Xtemp = Xscale(:,:,i);
    Xtemp(1:1+size(Xtemp,1):end) = 0; % diagonals are zero
    Xscale(:,:,i) = Xtemp;
end 
Xs = Xscale;
T = 298.15;
filled_ind = find(~isnan(Xs));
fillmethods = ["dia", "avg", "avr", "avc","uni"];

% vars for missing_parafac3
maxiter = 20000;
scale = 1;
center = 1;
conv = 1e-10;
winsorized_mse = 0;
orth = 1;
whichX = 'scale';

% loops, outer to inner: fill method, concentration, fn 
for iter = 5
    %find fill method to be used
    fillmethod = fillmethods(iter);

    %intialise the metrics 
    msefill= zeros(size(fns,2),size(concentrations,2));
    aardfill = zeros(size(fns,2),size(concentrations,2));
    wmsefill = zeros(size(fns,2),size(concentrations,2));

    %time the method 
    tic
    %j=0; % counter for intervals
    plotcount = 0;
    for c = [2,4]%:2:length(concentrations)
        conc = concentrations(c); 
        disp(conc)
        if strcmp(whichX, 'scale') 
            Xs = Xscale(:,:,c);
        else 
            Xs = Xsign(:,:,c);
        end 
        Xs = remove_nan2(Xs);
        dim = size(Xs);
        dim1 = dim(1);
        dim2 = dim(2);
        missing_ind = find(isnan(Xs));
        [row,col] = find(~isnan(Xs));
        filled_ind = find(~isnan(Xs));

            parfor i=1:length(fns)
                %i=i+1;
                fn = fns(i);
                disp('rank')
                disp(fn)
                %[U,D,V,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
                %[U,D,V,St,X_pred, iters]=missing_svd(Xs,fn,1,1,1e-8,1000);
                [X_pred,iters,F,err] = missing_parafac3(Xs,fn,maxiter,conv,scale,center, fillmethod, orth,mixtures, conc, whichX, T);
                %X_pred is the model predictions, not only missing values are
                %filled 
                % the actual model, not just filled values - can only use errors of filled values to find model error
                Xm= X_pred;
                msefill(i,c) = (sum((Xm(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
                wmsefill(i,c) = find_wmse(Xs(filled_ind), Xm(filled_ind), length(filled_ind));
                abserrorfill = abs(X_pred(filled_ind)-Xs(filled_ind));
                Xtemp = Xs(filled_ind);
                indices_use = find(Xtemp~=0); 
                aardfill(i,c) = sum(abserrorfill(indices_use)./Xtemp(indices_use))/length(Xtemp(indices_use));
            end %END FNS    
    end %END CONC 
    % Export results for this filling method 
    filenamesave = strcat('2wayPARAFAC-X',whichX,'-maxiter=10000-T=',num2str(T), '-fillmethod=', fillmethod ,'-orth=', num2str(orth),'-',date, '.3.mat');
    save(filenamesave)
end % END FILL METHODS  
toc

% %% Correlations
% % Change filename based on which corrs you want to find 
% fillmethods = ["dia", "avg", "avr", "avc","uni"];
% for iter =1:5
%     fillmethod = fillmethods(iter);
%     filename = strcat('2wayrunPARAFACPolySmall-0orth-298.15-',fillmethod,'-',date, '.mat');
%     load(filename)
%     disp(filename)
% 
%     corrFac = cell(length(conc_interval),N,2);
%     pvalFac = cell(length(conc_interval),N,2);
%     for r = 3:N
%         for j = 1:9
%         for i =1:2
%             [corrFac{j,r,i}, pvalFac{j,r,i}] = corr(Factors{j,r}{1,i});
%         end 
%         end 
%     end 
%     save(filename)
% end 

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
%%
Xf = filldata3(Xs,'uni',mixtures,conc,'sign');
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

