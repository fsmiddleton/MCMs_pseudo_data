%% 3-way array completion using PARAFAC and LOOCV 

%FS Middleton 2022/06/13

%% Import data 
clc
clear
T = 298.15; %K 
filename = strcat('TestSMALLHEMatrix13June',num2str(T),'.xlsx');
table = table2array(readtable(filename, 'Sheet', '0.1'));
conc_interval = 0.1:0.1:0.5;
dim1 = size(table,1);
dim2 = size(table,2);
dim3 = length(conc_interval);
X = nan(dim1, dim2, dim3);

for i = 1:length(conc_interval)
    table = table2array(readtable(filename, 'Sheet',num2str(conc_interval(i))));
    X(:,:,i) = table;
end 
% Analysis of the 4-way array 
percmiss = length(find(isnan(X)))/(dim1*dim2*dim3)*100;
percobs = length(find(~isnan(X)))/(dim1*dim2*dim3)*100;
%% Load data

interval = 0.05;
T = 298.15;
load(strcat('HEDataPolySmall-0.05-298.15.mat'))

%% Plt the missing data structure 
% these must have the same sizes as x
v=X;

xslice = 1:1:20;    % location of y-z planes
yslice = 1:2:20;     % location of x-z plane
zslice = 1:1:5;         % location of x-y planes
clf
slice(v,xslice,yslice,zslice)
xlabel('Component 1')
ylabel('Component 2')
zlabel('Composition')
% hm = HeatMap(HE_matrix);
% addXLabel(hm,'Component 1','FontSize',12);
% addYLabel(hm,'Component 2','FontSize',12);
% addZLabel(hm, 'Composition', 'FontSize',12);
% view(hm)
% histogram(X)
% xlabel('Value in the matrix X')
% ylabel('Frequency')

%% Find the correct number of factors 
prefixfilename = ('HEData3wayPolySmall');
Temps = 298.15; %k
interval = [0.01; 0.02; 0.04; 0.05];

tic
for iter = 4
    intervalloop = interval(iter,1);
    % Import array from .mat file
    filename = strcat(prefixfilename, '-',num2str(intervalloop),'-', num2str(Temps), '.mat');
    load(filename)
    conc_interval = intervalloop:intervalloop:(0.5-intervalloop);
    X = HE_data_sparse(:,:,1:length(conc_interval));
    %define things necessary in the completion 
    N=12; % used in a loop therefore can be a vector 
    max_iter = 10000;
    conv = 1e-10;
    scale = 1;
    center = 1;
    fillmethod = 'avg';

    %Find the correct number of factors 
    
    %  number of factors maximum to try 
    

    filled_ind = find(~isnan(X));
    % for each % missing 
    minmse = zeros(1,length(Temps));
    minaapd = zeros(1,length(Temps)); % average absolute percent deviation = average(abs(Xpred-Xtrue)/Xtrue))
    numberOfFactors = zeros(1,length(Temps));
    dof = zeros(1,length(Temps));
    mse = zeros(N,length(Temps));
    msefill = zeros(N,length(Temps));
    RADfill = zeros(N,length(Temps));
    aapdfill = zeros(N,length(Temps));
    indexsortmse=zeros(N,length(Temps));
    c = zeros(N,length(Temps));
    coreconsistency = zeros(N,1);
    randomstart = zeros(N,length(Temps));

    %metrics defined
    smse = zeros(N,1);
    fit = zeros(N,1);
    it = zeros(N,1);
    error = zeros(N, length(filled_ind));%missing_ind is biggest for 90% missing data, as initialisation was
    errorfill = zeros(N,length(filled_ind));

    % define other needed vars  
    dim = size(X);
    X_model = zeros(dim(1),dim(2),dim(3),N);
    count = 0; % index for the missing data arrays
    Factors = cell(N,1);
    % time the process 


    for temperature = Temps
        disp('Temperature')
        disp(temperature)
        count = count+1;

        no_filled = temperature/100*dim(1)*dim(2)*dim(3);
        %remove nan slabs from the 3-way array 
        %[X,dim, znan, ynan, xnan]=remove_nan3(X);
        % Find the correct number of factors for this percent missing 
        for n =N % factors (== principal components) 
            % initialise random number generator 
            randind = 1;
            rng(randind, 'twister')

            disp('n')
            disp(n) 
            %perform INDAFAC with missing values on the matrix with one entry missing
            %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
             [X_pred,iter,F,err] = missing_parafac3(X,fn,max_iter,conv,scale,center,  fillmethod);
            Xm = X_pred;
            %mixture
            Factors{n,1}=F;
            errorfill(n,:) = (X(filled_ind)-Xm(filled_ind))'; % actual prediction error
            X_model(:,:,:,n) = X_pred;
            % averages of metrics 
            msefill(n,count) = sum(errorfill(n,:).^2)/length(filled_ind);

           
            RAD = errorfill(n,:)./(X(filled_ind))';
            RAD = RAD(find(~isinf(RAD)));
            RADfill(n,count)= sqrt((sum(RAD).^2)/length(filled_ind));
            randomstart(n,count) = randind;
            abserrorfill = abs(errorfill(n,:));
            %built in metric
            %fit(n,count) = D.fit; 
            error(n,:) = (X(filled_ind)-Xm(filled_ind))';
            mse(n,count) = sum(error(n,:).^2)/length(filled_ind);
            %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
            %c(n,count)= corcond(Xm,F,[],1);
        end
        % Find true number of factors for this % missing 
        %Find the optimal rank prediction
        [sortmse, indexsortmse(:,count)]=sort(mse(:,count)); % sorts in ascending order 
        m = 1; % counter for finding consisent number of factors with lowest mse 
        minmse(count) = sortmse(m);
        numberOfFactors(count) = indexsortmse(m,count);
        cmin = coreconsistency(numberOfFactors(count));
        dof(count) = dim(1)*dim(2)*dim(3)-numberOfFactors(count)*(dim(1)+dim(2)+dim(3)-2);
        minaapd(count) = aapdfill(numberOfFactors(count),count);
        corrFac = cell(N,3);
        pvalFac = cell(N,3);
        filenametemp = strcat('3wayrunINDAFACPolyCenter+Scale-',fillmethod,'-int-', num2str(intervalloop), '-0orth-',num2str(T),'-',date, '.mat');
        save(filenametemp)
    end
end 
toc % end timer 



%% LOOCV for the correct number of factors 
% load data 
prefixfilename = ('HEData3wayPolySmall');
Temps = 298.15; %k
interval = 0.05;
filename = strcat(prefixfilename, '-',num2str(interval),'-', num2str(Temps), '.mat');
load(filename)
conc_interval = interval:interval:(0.5-interval);
X = HE_data_sparse(:,:,1:length(conc_interval));
% everything necessary for this to run 
rng(2,'twister')
%Choose number of factors
numberOfFactorsLOOCV = 6:9;% = numberOfFactors(missing==missingLOOCV);
fillmethod = 'dia';
T = 298.15;
max_iter = 3000;
conv = 1e-10;
center = 1;
scale = 1;
% find filled and missing indices
[row,col] = find(~isnan(X(:,:,1))); % one 2-way slice 
filled_ind= find(~isnan(X(:,:,1))); 
missing_ind = find(isnan(X)); 
[i,j,k]=findfill3(X);
dim = size(X); 


% define other needed vars. These are rewritten on each iteration of number
% of factors in a component 
N = length(filled_ind);
smseN = zeros(1,N);
fitN = zeros(1,N);
itN = zeros(1,N);
errorN = zeros(length(conc_interval),N);
RAD = zeros(length(conc_interval),N);
mseN = zeros(1,N);
cN = zeros(1,N);
coreconsistencyN = zeros(1,N);
X_predN = zeros(dim(1),dim(2),dim(3),N);
n = 1; %just a counter 
% LOOCV

% time the LOOCV
tic

for factors=numberOfFactorsLOOCV
    errorN = zeros(length(conc_interval),N);
    RAD = zeros(length(conc_interval),N);
    Factors = cell(7,N);
    cN = zeros(1,N);
    X_point = zeros(length(filled_ind),1);
    parfor count = 1:length(filled_ind) %filled_linear_ind must be a row vector  
        % temporary X
        X_b = X;
        %count=count+1;
        disp(count)
        X_b(row(count),col(count),:) = nan;% remove a point from Xs
        if   find(~isnan(X_b(i(count),:,k(count))))&  find(~isnan(X_b(i(count),:,k(count)))) & (X_b(row(count),col(count),1))~=0%ensure at least one value in each column and row
            %count=count+1; % allowed, therefore the counter is increased 
            %perform INDAFAC with missing values on the matrix with one entry missing
            %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
            [X_pred,iter,F,err] = missing_parafac3(X_b,factors,max_iter,conv,scale,center, fillmethod);
            Xm = X_pred;
            Factors{factors, count} = F;
            %predicted data point 
            X_point(count) = Xm(count);
            
            errorN(:,count) = X(row(count),col(count),:)-Xm(row(count),col(count),:);
            RAD(:,count) = errorN(:,count)./reshape(X(row(count),col(count),:),size(errorN(:,count)));
            tempvalue = 0;
            %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
            cN(n,count)= corcond(Xm,F,[],1);


        end 
    end %end LOOCV 
%     for r = factors
%         for i =1:3
%             [corrFac{r,i}, pvalFac{r,i}] = corr(Factors{r,1}{1,i});
%         end 
%     end
    
    % averages of metrics for LOOCV 
    msebest = sum(errorN(n,:).^2)./length(errorN(n,:));
    RAD=RAD(find(~isinf(RAD)));
    RAD_overall_fill = sum(RAD(n,:).^2)./length(RAD(n,:));
    coreconsistencybest = sqrt(sum(cN(n,:).^2)./length(cN(n,:)));
    X_pred_LOOCV=X;
    X_pred_LOOCV(filled_ind)=X_point;
    filenametemp = strcat('3wayrunPARAFACPolySmall-LOOCV-N=',num2str(factors), '-',num2str(Temps),'-interval-', num2str(interval), '-', fillmethod,date, '.mat');
    save(filenametemp)
end 
toc 

%% LOOCV old 
% must be run after the correct nnumber of factors is found 
rng(2,'twister')
%Choose number of factors
numberOfFactorsLOOCV = 1;% = numberOfFactors(missing==missingLOOCV);
fn = numberOfFactorsLOOCV;
[row,col] = find(~isnan(X(:,:,1))); % one 2-way slice 
filled_ind= find(~isnan(X(:,:,1))); 
missing_ind = find(isnan(X));
% find filled indices 
[i,j,k]=findfill3(X);
dim = size(X);
% define metrics here that need to be used 
N = length(filled_ind);
smseN = zeros(1,N);
fitN = zeros(1,N);
itN = zeros(1,N);
errorN = zeros(length(conc_interval),N);
AAD = zeros(length(conc_interval),N);
mseN = zeros(1,N);
cN = zeros(1,N);
coreconsistencyN = zeros(1,N);
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%
% define how to use method again 
center = 1;
scale = 1;
maxiter = 1000;
modeINDAFAC =1;
method = 'Broindafac';

% define other needed vars 
X_predN = zeros(dim(1),dim(2),dim(3),N);
n = 1; % factors (== principal components) 
% LOOCV
count=0;% counter for LOOCV
ind = 0; % counter for filled indices 

% time the LOOCV
tic
for filled_index = filled_ind' %filled_linear_ind must be a row vector  
    % temporary X
    X_b = X;
    ind=ind+1;
    disp(ind)
    X_b(row(ind),col(ind),:) = nan;% remove a point from Xs
    if   find(~isnan(X_b(i(ind),:,k(ind))))&  find(~isnan(X_b(i(ind),:,k(ind)))) & (X_b(row(ind),col(ind),1))~=0%ensure at least one value in each column and row
        count=count+1; % allowed, therefore the counter is increased 
        %perform INDAFAC with missing values on the matrix with one entry missing
        %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
        [X_pred,iter,F,err] = missing_parafac3(X_b,fn,maxiter,conv,scale,center, fillmethod);
        Xm = X_pred;
        %built in metric
        %fitN(n,count) = D.fit;
        %own metrics 
        % fix this 
        errorN(:,count) = X(row(ind),col(ind),:)-Xm(row(ind),col(ind),:);
        AAD(:,count) = errorN(:,count)./reshape(X(row(ind),col(ind),:),size(errorN(:,count)));
        %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
        cN(n,count)= corcond(Xm,F,[],1);
    end 
end %end LOOCV 
% averages of metrics for LOOCV 
msebest = sum(errorN(n,:).^2)./length(errorN(n,:));
AAD=AAD(find(~isinf(AAD)));
AAD_overall = sum(AAD(n,:).^2)./length(AAD(n,:));
coreconsistencybest = sqrt(sum(cN(n,:).^2)./length(cN(n,:)));

toc 


%% 
function [X,dim, znan, ynan, xnan]=remove_nan3(X)
    % saves the columns and rows with only Nan values and removes them from
    % the 3-way array 
    % Input 
    % X = 3-way array 
    % Output
    % X = 3-way array without slabs of only nan values 
    % dim = dimensions of X outputted 
    % znan = index of z coordinates that contained only nan values 
    % ynan = index of y coordinates that contained only nan values 
    % xnan = index of x coordinates that contained only nan values 
    dim = size(X);
    % z slabs
    t1 = all(isnan(X),[1,2]);
    t1 = reshape(t1,[1,dim(3)]); %logical 
    r1=find(t1);
    X(:,:,r1)=[];

    % x slabs 
    t2 = all(isnan(X),[2,3]);
    t2 = reshape(t2,[1,dim(1)]); %logical 
    r2=find(t2);
    X(:,:,r2)=[];

    % y slabs 
    t3 = all(isnan(X),[1,3]);
    t3 = reshape(t3,[1,dim(2)]); %logical 
    r3=find(t3);
    X(:,:,r3)=[];

     %new X 
    dim = size(X);
    znan=r1;
    xnan=r2;
    ynan=r3;
end 

function [X, Xtrue, Factors] = importX(filename, dim)
    Factors = readtable(filename, 'Sheet', 'Factors');
    X=zeros(dim(1),dim(2),dim(3));
    Xtrue=zeros(dim(1),dim(2),dim(3));
    for i =1:dim(1)
        Ttrue = readtable('ToyProblemData3DFull.xlsx','Sheet', num2str(i));
        Xtrue(i,:,:)= table2array(Ttrue);
        T = readtable(filename, 'Sheet', num2str(i));
        X(i,:,:)=table2array(T);
    end
end

function [F,D, X_pred]=missing_indafac(X,fn,mode, center,scale,conv,max_iter,method)
% Input 
% X = data array 
% fn = number of factors 
% mode = mode of unfolding used 
% center  =1 center data 
% scale  =1 scale data 
%         =0 do not scale data  
% conv = convergence criterion = difference between new and old 
% method ='parafac' or 'indafac' or 'Broindafac'
 
% Output 
% F = factors (A,B,C)
% D = diagnostics 
%     D.it = iterations 
%     D.fit = fit
%     D.error = error
%     D.cor = core consistency 
%     D.stop = reason for stopping (max iterations or conv reached)
% X_pred = predicted array 

% INDAFAC is used to find the factors of a 3-way array for a given number
% of factors (fn), with options as to how to unfold the matrix for centering and scaling, 
% as well as whether to center and scale and when to stop iterations 


    dim=size(X);
    [i,j,k]=findnan3(X); % missing indices 
    if method == 'Broindafac'
        % use their indafac algorithm, not my method 
        Fi = ini(X, fn, 1); % initialise the three loadings 
        [F,D]=parafac(X,fn);
        Xm = nmodel(F);
        X_pred = X;
        % fill misssing values 
        for d = 1:length(i)
            X_pred(i(d),j(d),k(d))= Xm(i(d),j(d),k(d));
        end
    else 
        % this only uses indafac to find the loadings, with my own
        % centering and updating 
        % much longer to run algorithm 
        if length(i)>1 % there is missing data 
            Xf = filldata3(X); % fill slices with averages  
            SS = sumsquare3(Xf,i,j,k); % sum of squares of missing values 
            X_pred = X;
            f=2*conv;
            iter = 1;
            while iter<max_iter && f>conv
                SSold = SS;
                % preprocess = scale and then center 

                if scale ==1 || center ==1
                    [x1,x2,x3]=nshape(Xf);
                    %unfold
                    if mode ==1
                        Xp=x1;
                    elseif mode ==2
                        Xp = x2;
                    else 
                        Xp = x3;
                    end
                    %scale 
                    if scale ==1 
                        sj = sqrt(sum(sum((Xp).^2)));
                        Xp = Xp/sj;
                    end 
                    % center
                    if center ==1
                        mx = mean(Xp);
                        % Center X across columns
                        Xp = Xp-mx*ones(size(mx,2),1);
                    end 
                    % reshapes back to 3-way array 
                    Xc=reshape(Xp,size(X));
                else
                    Xc = Xf;
                end 

                % Find factors 
                % initial guess for factors 

                Fi = ini(Xc, fn, 1); % initialise the three loadings
                alphabet = 'ABCD'; 

                for alph=1:3 %each loading for each dimension of the data
                    eval([alphabet(alph) ,'= Fi{alph};']);
                end 
                % indafac or parafac step 
                if strcmp(method,'indafac') 
                    [F,D]=INDAFAC(Xc,fn,Fi);
                    Xm = nmodel(F);
                else 
                    [F, ~]=parafac(Xc,fn);
                    Xm = nmodel(F);
                end 

                % post process = uncenter, unscale 
                if scale ==1 || center ==1
                    [x1,x2,x3]=nshape(Xm);
                    %unfold
                    if mode ==1
                        Xp=x1;
                    elseif mode ==2
                        Xp = x2;
                    else 
                        Xp = x3;
                    end
                    % uncenter
                    if center ==1
                        Xp = Xp+mx*ones(size(mx,2),1);
                    end 
                    % unscale 
                    if scale ==1 
                        Xp = Xp*sj;
                    end 
                    % reshapes back to 3-way array 
                    X_pred=reshape(Xp,size(X));
                else
                    X_pred = Xm;
                end  

                % fill misssing values 
                for d = 1:length(i)
                    Xf(i(d),j(d),k(d))= X_pred(i(d),j(d),k(d));
                end 
                % check for convergence 
                SS = sumsquare3(Xf,i,j,k); % sum of squares of missing values 
                f = abs(SS-SSold)/(SSold);
                iter = iter+1;
            end
            X_pred = Xf;
        else 
            % no missing data
            disp('No missing data')
        end %end if  else
    end
end % end function 

function [X_filled, missing] = filldata3(X, method,mixtures,concintervalarray)
    % filldata3 fills a 3-way array with the arithmetic mean of the other
    % values in the array
    % Input 
    % X = data array 
    % method = 'avg' - column and row averages 
    %       or 'uni' - unifac predictions 
    % mixtures = components ordered according to their linear place on the 2-way array 
    % concinterval = concentration interval over which to fill data 
    % Output 
    % X_filled = filled data array
    % missing = linear indices of missing data 

    % fill as done for PCA with SVD (averages of each dimension)
    [i, j, k]=findnan3(X);% returns rows and columns with nan elements
    dim = size(X);
    missing = isnan(X); %linear indices 

    X_filled = X;
    % fill missing values with zeros 
    for ind = 1:length(i)
        X_filled(i(ind),j(ind), k(ind))=0;
    end 
    
    if strcmp(method,'avg')
    % find all means 
        mean1 = sum(X_filled,1)./(ones(1,dim(2),dim(3))*dim(1)-sum(missing,1)); % (1,j,k) dim1 
        mean2= sum(X_filled,2)./(ones(dim(1),1,dim(3))*dim(2)-sum(missing,2)); % (i,1,k)  dim2 
        mean3 = sum(X_filled,3)./(ones(dim(1),dim(2),1)*dim(3)-sum(missing,3));% (i,j,1) dim3 
        %replace nan means with 0 
        mean1(find(isnan(mean1)))=0;
        mean2(find(isnan(mean2)))=0;
        mean3(find(isnan(mean3)))=0; 
        for ind =1:length(i)
           % for all NaN elements that exist, loop through them to replace
            X_filled(i(ind),j(ind), k(ind))=(mean1(1,j(ind), k(ind))+mean2(i(ind),1, k(ind))+mean3(i(ind),j(ind),1))/3;
        end 
    elseif strcmp(method,'dia')
        X_filled = X;
        for d = 1:dim(3)
            Xtemp = X(:,:,d);
            [m,n]=size(Xtemp);
            missing_ind = (isnan(Xtemp));
            [i, j]=find(isnan(Xtemp));% returns rows and columns with nonzero elements
            X_filledtemp=Xtemp;
            X_filledtemp(missing_ind)=0; %fill NaN values with 0   
            % for all NaN elements that exist, loop through them to replace with means 
            
            % lower and upper triangular parts of the array
            LX = tril(X_filledtemp,-1);
            UX = triu(X_filledtemp, 1);
            % averages for lower and upper
            mean_rowL = sum(LX,2)./(ones(1,m)*n-sum(missing_ind,2)/2);
            mean_colL = sum(LX,1)./(ones(1,n)*m-sum(missing_ind,1)/2);
            mean_rowU = sum(UX,2)./(ones(1,m)*n-sum(missing_ind,2)/2);
            mean_colU = sum(UX,1)./(ones(1,n)*m-sum(missing_ind,1)/2);
            
            for ind =1:length(i) 
                if i(ind)>j(ind) % lower diagonal
                    X_filledtemp(i(ind),j(ind))=(mean_colL(j(ind))+mean_rowL(i(ind)))/2;
                else 
                    X_filledtemp(i(ind),j(ind))=(mean_colU(j(ind))+mean_rowU(i(ind)))/2;
                end 
            end
            X_filled(:,:,d)=X_filledtemp;
        end 
    else %method ='uni'
        % linearise X to a column vector - entries match mixtures
        X1 = reshape(X(:,:,1),[dim(1)*dim(2),1]); % column vector used for checks and indices
        Xtemp = reshape(X,[dim(1)*dim(2),dim(3)]); % reshaped to column vectors with each column containing a concentration interval
        
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
        X_filled = reshape(Xtemp,[dim(1), dim(2), dim(3)]);
        %fills remaining few Nan values with averages 
        X_filled = filldata3(X_filled,'avg',mixtures,conc_interval);
    end 
end 

function [i,j,k]=findfill3(X)
% Input 
% X = data array 
% Output 
% i, j, k = indexes of existing values in X

% findnan3 finds the nan indices of a 3-way array 
    dim3 = size(X,3);
    i=[];
    j=[];
    k=[];
    for d = 1:dim3
        % done per z slice 
        Xtemp = X(:,:,d);
        [itemp, jtemp]= find(~isnan(Xtemp));
        i = [i; itemp];
        j = [j;jtemp];
        ktemp = ones(length(itemp),1)*d;
        k = [k;ktemp];
    end 
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

function [SS]=sumsquare3(X,i,j,k)
% Input 
% X = data array 
% [i,j,k] indices of values of interest
% Output 
% SS = sum of the squares of the values of interest 
    tempX=[];
    for d = 1:length(i)
        tempX=[tempX;X(i(d),j(d),k(d))];
    end 
    SS = sum(tempX.^2);
end 


