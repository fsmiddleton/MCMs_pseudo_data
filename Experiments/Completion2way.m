%% 2-way arrays 
% Francesca Middleton, 2022-03-02

%% Import toy data 

clc
clear

import = 1; 
% import = 0 to create a low rank matrix, 
% import = 1 to fetch an already created low rank matrix (2-way) 
% import = 2 to import excess enthalpy data (2-way),

filename = 'RandomMatrixNoise.xlsx';
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
    T1 = readtable('TestSMALLHEMatrix13June298.15.xlsx', 'Sheet', '0.1'); 
    X = table2array(T1);%only use T1, small enough to manage
    
else
    disp('Please choose 0, 1, 2 to use data')
end 
%% Import or create data matrix for testing 3-way 
clc
clear
% Import excess enthalpy data for one matrix 
filename = 'TestHEMatrixSMALL15June298.15.xlsx';
T1 = readtable(filename, 'Sheet', '0.1'); 
X = table2array(T1);%only use T1, small enough to manage
dim = size(X);
dim1=dim(1);
dim2=dim(2);
conc_interval = 0.1:0.1:0.9;
percmiss = length(find(isnan(X)))/(dim1*dim2)*100;
percobs = length(find(~isnan(X)))/(dim1*dim2)*100;
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
%create the 3-way array with filled values - to be slices 
dim1 = size(X,1);
dim2 = size(X,2);
dim3 = length(conc_interval);
X3 = nan(dim1, dim2, dim3);

for i = 1:length(conc_interval)
    table = table2array(readtable(filename, 'Sheet',num2str(conc_interval(i))));
    X3(:,:,i) = table;
end 

Xf_ini = filldata3(X3,'uni',mixtures,conc_interval); 

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
filename = 'TestSMALLHEMatrix13June298.15.xlsx';
fn = min([dim1,dim2]);
y1 = zeros(fn,length(concentrations));
y2 = zeros(fn,length(concentrations));

%necessary for plots
indexplot=20;

i=0;
for c = concentrations
    i=i+1;
    %Import data 
    T = readtable(filename, 'Sheet', num2str(c));
    Xs=table2array(T);
    %SVD
    [~,D,~,~,~, ~]=missing_svd(Xs,fn,1,1,1e-8,1,0);
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
%% LOOCV to find rank using mse  
%time the code 
tic
% declare missing % used and the file from which to extract information 
fns =7:1:10;
concentrations=0.1:0.1:0.5;
Xs=X;
missing_ind = find(isnan(Xs));
% largest length of filled_linear_ind is for 0.3
filled_ind = find(~isnan(Xs));

% vars for parafac
maxiter = 10000;
scale = 1;
center = 1;
conv = 1e-8;
winsorized_mse = 0;
fillmethod = 'avg';

% declare vars for analysis 
mse_LOOCV = zeros(length(concentrations),length(fns));
wmse_LOOCV = zeros(length(concentrations),length(fns));
RAD_LOOCV = zeros(length(concentrations),length(fns)); % relative absolute deviation 
min_mse_LOOCV = zeros(length(concentrations),1);
min_wmse_LOOCV = zeros(length(concentrations),1);
X_pred_LOOCV = zeros(dim1,dim2,length(filled_ind));
LOOCV_removed_col = zeros(length(filled_ind),1);
LOOCV_removed_row = zeros(length(filled_ind),1);
SVs = zeros(dim1,length(concentrations));

% Error bars - LOOCV 
% import the relevant sparse matrix from the spreadsheet
ind = 0;

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
    disp('Concentration of compound 1')
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
                
                [X_pred,iters,F,err] = missing_parafac3(X_b,fn,maxiter,conv,scale,center,fillmethod);
                X_pred_LOOCV(:,:,k)=X_pred;
                error_LOOCV(fn,k) = Xs(filled_index)-X_pred(filled_index);
                
                Xm_boot(c, fnind, k) = X_pred(filled_index);
                if Xs(filled_index)~=0
                    RAD(fn,k) = error_LOOCV(fn,k)/Xs(filled_index);
                end
                
            end
        end
        % mse for this composition and rank 
        mse_LOOCV(c,fnind)= sum(error_LOOCV(fn,:).^2)/length(error_LOOCV(fn,:));
        wmse_LOOCV(c, fnind) = find_wmse_error(error_LOOCV(fn,:), length(filled_ind'));
        %absolute average deviation
        RAD_LOOCV(c,fnind) = (sum(abs(RAD(fn,:))))/length(RAD(fn,:));
    end % END FN
    filenamesave = strcat('2wayPARAFACSmall-LOOCV-maxiter=10000-c=', num2str(c),  date, '.mat');
    save(filenamesave)
    % find the optimal rank 
    if winsorized_mse ==1
        fn_LOOCV(c) = find(wmse_LOOCV(c,:) == min(wmse_LOOCV(c,:)));
    else 
        fn_LOOCV(c) = find(mse_LOOCV(c,:) == min(mse_LOOCV(c,:)));
    end
    
    min_mse_LOOCV(c) = mse_LOOCV(c,fn_LOOCV(c));
    min_wmse_LOOCV(c) = wmse_LOOCV(c,fn_LOOCV(c));
    % use the best rank to find error bars 
end 
toc 

%% Plot errors 
clf
conc=0.1;
c=conc*10;
plot(fns, (RAD_LOOCV(c,:)))

%% Find optimal rank of the model for the filled data 

% Import other data for most missing to allow variables to be declared
c=30;
T = readtable(filename, 'Sheet', num2str(c/100));
Xs=table2array(T);

%necessary vars 
dim = size(Xs);
dim1=dim(1);
dim2=dim(2);
missing_ind = find(isnan(Xs));

%vars to loop through 
fillmethods = ["dia", "avg", "avr", "avc","uni"];
fns =3:1:12;
concentrations = 50:10:90;

%vars governing the parafac algorithm 
maxiter =10000;
conv = 1e-8;
scale = 1;
center = 1;
reorder = 0;
winsorized_mse = 1;
fillmethod = 'avg';

% loops: fill method, concentration, fn 
for iter = 1:2
    
    % can also get Rho, Lambda, EV% and Max Gr
    % EV must be close to 100% - explained variation 
    dim = size(X);
    %Find the correct number of factors 
    Temps = 298.15;%K
    %  number of factors maximum to try 
    
    fillmethod = fillmethods(iter);
    disp(fillmethod)

    %intialise the metrics to analyse each sparsity and its final rank found 
    % vectors, one for each matrix that is sliced 
    minmse_miss = zeros(size(concentrations,2),1); 
    minwmse_miss = zeros(size(concentrations,2),1);
    minmse_fill = zeros(size(concentrations,2),1);
    minaapd = zeros(length(concentrations),1);
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
    aapdfill = zeros(size(fns,2),size(concentrations,2));
    wmse_miss = zeros(size(fns,2),size(concentrations,2));
    wmsefill = zeros(size(fns,2),size(concentrations,2));
    R = zeros(size(fns,2),size(concentrations,2));
    cumulative_contribution = zeros(length(concentrations),length(fns));
    SVs = zeros(length(fns),length(concentrations));

    %For gavish 
    ranks_attempted = zeros(1000,length(concentrations));
    %time the method 
    tic
    j=0; % counter for intervals
    plotcount = 0;
    for c = concentrations

        j=j+1; % for sparsity
        disp('missing')
        disp(c)
        %Import sparse matrix for toy problem 
        T = readtable(filename, 'Sheet', num2str(c/100));
        Xs = table2array(T);
        Xs = remove_nan2(Xs);
        dim = size(Xs);
        dim1 = dim(1);
        dim2 = dim(2);
        missing_ind = find(isnan(Xs));
        [row,col] = find(~isnan(Xs));
        filled_ind = find(~isnan(Xs));

        % Find the optimal rank for this matrix

        
            % Iterative PCA with wMSE or MSE used to find rank

            %Find rank by minimising the mse or wmse 
            i=0;
            for fn=fns
                i=i+1;
                disp('rank')
                disp(fn)
                %[U,D,V,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
                %[U,D,V,St,X_pred, iters]=missing_svd(Xs,fn,1,1,1e-8,1000);
                [X_pred,iters,F,err] = missing_parafac3(Xs,fn,maxiter,conv,scale,center, fillmethod);
                %X_pred is the model predictions, not only missing values are
                %filled 
                % the actual model, not just filled values - can only use errors of filled values to find model error
                Xm= X_pred;
                msefill(i,j) = (sum((Xm(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
                wmsefill(i,j) = find_wmse(Xs(filled_ind), Xm(filled_ind), length(filled_ind));
                abserrorfill = abs(X_pred(filled_ind)-Xs(filled_ind));
                aapdfill(i,j) = sum(abserrorfill./X(filled_ind))/length(filled_ind);
    %             Cyt = corrcoef(X_pred(missing_ind),Xs(missing_ind));
    %             R(i,j)=sqrt(Cyt(2,1));
            end %END FNS

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
            minaapd(j) = aapdfill(min_index,j);
 
    end %END CONC 
    % Export results 
    filenamesave = strcat('2wayPARAFACPolySmall-2-fill=', fillmethod,'-', date, '.mat');
    save(filenamesave)
end % END FILL METHODS  
toc

%% Correlations
% Change filename based on which corrs you want to find 
fillmethods = ["dia", "avg", "avr", "avc","uni"];
for iter =1:5
    fillmethod = fillmethods(iter);
    filename = strcat('2wayrunPARAFACPolySmall-0orth-298.15-',fillmethod,'-',date, '.mat');
    load(filename)
    disp(filename)

    corrFac = cell(length(conc_interval),N,2);
    pvalFac = cell(length(conc_interval),N,2);
    for r = 3:N
        for j = 1:9
        for i =1:2
            [corrFac{j,r,i}, pvalFac{j,r,i}] = corr(Factors{j,r}{1,i});
        end 
        end 
    end 
        
    save(filename)
end 

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
function [X_filled, missing] = filldata3(X, method,mixtures,concintervalarray)
    % filldata3 fills a 3-way array with the arithmetic mean of the other
    % values in the array. Used for a 2-way array here
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
    % returns rows and columns with nan elements
    dim = size(X);
    
    if (size(dim,2))==2
        dim(3)=1;
    end 
    missing = isnan(X); %linear indices 

    X_filled = X;
     
    
    if strcmp(method,'avg')
    % f
        % in 2-wy array, do 2-way avg fills
        X_filled = zeros(dim);
        for i = 1:dim(3)
            X_filled(:,:,i) = fill_data(X(:,:,i)); 
        end 
    elseif strcmp(method, 'avc')
        % column averages are used only 
        [m,n]=size(X);
        missing_ind = (isnan(X));
        [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
        X_filled=X;
        X_filled(missing_ind)=0; %fill NaN values with 0
        mean_col = sum(X_filled,1)./(ones(1,n)*m-sum(missing_ind,1));   
        % for all NaN elements that exist, loop through them to replace with means 
        for k =1:length(i) 
            X_filled(i(k),j(k))=mean_col(j(k));
        end  
    elseif strcmp(method, 'avr')
    % row averages are used only 
        [m,n]=size(X);
        missing_ind = (isnan(X));
        [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
        X_filled=X;
        X_filled(missing_ind)=0; %fill NaN values with 0
        mean_row = sum(X_filled,2)./(ones(1,m)*n-sum(missing_ind,2));   
        % for all NaN elements that exist, loop through them to replace with means 
        for k =1:length(i) 
            X_filled(i(k),j(k))=mean_row(i(k));
        end
    elseif strcmp(method, 'dia') % lower and upper diagonals are used in their averages 
        [m,n]=size(X);
        missing_ind = (isnan(X));
        [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
        X_filled=X;
        X_filled(missing_ind)=0; %fill NaN values with 0   
        % for all NaN elements that exist, loop through them to replace with means 
        
        % lower and upper triangular parts of the array
        LX = tril(X_filled,-1);
        UX = triu(X_filled, 1);
        % averages for lower and upper
        mean_rowL = sum(LX,2)./(ones(1,m)*n-sum(missing_ind,2)/2);
        mean_colL = sum(LX,1)./(ones(1,n)*m-sum(missing_ind,1)/2);
        mean_rowU = sum(UX,2)./(ones(1,m)*n-sum(missing_ind,2)/2);
        mean_colU = sum(UX,1)./(ones(1,n)*m-sum(missing_ind,1)/2);
        
        for k =1:length(i) 
            if i(k)>j(k) % lower diagonal
                X_filled(i(k),j(k))=(mean_colL(j(k))+mean_rowL(i(k)))/2;
            else 
                X_filled(i(k),j(k))=(mean_colU(j(k))+mean_rowU(i(k)))/2;
            end 
        end
    else %method ='uni'
        % linearise X to a column vector - entries match mixtures
        X1 = reshape(X(:,:,1),[dim(1)*dim(2),1]); % column vector used for checks and indices
        Xtemp = reshape(X,[dim(1)*dim(2), dim(3)]); % reshaped to column vectors with each column containing a concentration interval
        
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
        
            
        concindices = find(strcmp(conc_array,conc_unifac));
        concindices2 = find(strcmp(conc_array2,conc_unifac));
        
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
        X_filled = reshape(Xtemp,[dim(1), dim(2), dim(3)]);
        %fills remaining few Nan values with averages 
        X_filled = filldata3(X_filled,'avg',mixtures,conc_interval);
    end 
end 

function  [fn,msefill,wmsefill,R2, X_pred,cutoff,SVs,ranks_iterations] = solveGavish(Xs, dim1, dim2, conv, iter)
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
    ranks_iterations = zeros(iter,1);
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
    
    if dim1<dim2
        fn = dim1;
    else 
        fn=dim2;
    end 
    
    it=0;
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
        ranks_iterations(it)=fn;
    end 
    disp('f')
    disp(f)
    disp('iter')
    disp(it)
    %metrics 
     % solve for the correct number of factors 
    [~,~,D,St,X_pred,~,~] = missing_svd(Xs,fn,1,1,1e-3,1);
    Xm = St*D';
    msefill = (sum((Xm(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
    wmsefill = find_wmse(Xs(filled_ind), Xm(filled_ind), length(filled_ind));
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

function [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd(X,fn,center,scale,conv,max_iter, use_missing)
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
            if scale ==1 
                sj = sqrt(sum(sum((Xfilled).^2)));
                Xfilled = Xfilled/sj;
            end 
            if center ==1
                mx = mean(Xfilled); %columnwise mean 
                Xc = Xfilled-ones(m,1)*mx; %centering of the data done - subtract mean from each entry in a column 
            else 
                Xc=Xfilled;
            end 
            %perform SVD
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
            SS = sum(sum((X_pred(indices)).^2));
   
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

function [wmse]=find_wmse_error(errors, count)
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





