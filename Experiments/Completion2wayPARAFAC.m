%% 2-way array completion using LOOCV 

%FS Middleton 2022/06/13

%% Import data 
clc
clear
T = 288.15; %K 
Temps = [298.15;243.15; 253.15; 263.15; 273.15; 283.15; 288.15; 290.15; 293.15; 296.15; 298.15; 303.15; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15; 363.15];

filename = strcat('HEMatrixPoly16June',num2str(T),'.xlsx');
table = table2array(readtable(filename, 'Sheet', '0.1'));
conc_analysis = 0.1;
dim1 = size(table,1);
dim2 = size(table,2);
%dim3 = length(conc_interval);
X = nan(dim1, dim2);


table = table2array(readtable(filename, 'Sheet',num2str(conc_analysis)));
X(:,:) = table;

% Analysis of the 2-way array 
percmiss = length(find(isnan(X)))/(dim1*dim2)*100;
percobs = length(find(~isnan(X)))/(dim1*dim2)*100;
%parpool
%% Find the best rank for each 

%assumes no missing values in any row/column
%Find the correct number of factors 
conc_interval = 0.1:0.1:0.9;%K
%  number of factors maximum to try 
N=12;
maxiter=10000;

missing_ind = find(isnan(X));
[row,col]=find(~isnan(X));
filled_ind = find(~isnan(X));
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%

% for each concentration and possible number of factors
minmse = zeros(1,length(conc_interval));
minaapd = zeros(1,length(conc_interval)); % average absolute percent deviation = average(abs(Xpred-Xtrue)/Xtrue))
numberOfFactors = zeros(1,length(conc_interval));
dof = zeros(1,length(conc_interval));
mse = zeros(N,length(conc_interval));
msefill = zeros(N,length(conc_interval));
RADfill = zeros(N,length(conc_interval));
aapdfill = zeros(N,length(conc_interval));
indexsortmse=zeros(N,length(conc_interval));
corecon = zeros(N,length(conc_interval));
coreconsistency = zeros(N,1);
randomstart = zeros(N,length(conc_interval), length(filled_ind));

%metrics defined
smse = zeros(N,1);
fit = zeros(N,1);
it = zeros(N,1);
error = zeros(N, length(conc_interval),length(filled_ind));%missing_ind is biggest for 90% missing data, as initialisation was
errorfill = zeros(N, length(conc_interval),length(filled_ind));
RAD = zeros(N, length(conc_interval),length(filled_ind));
% metrics to compare to truth 
msemiss = zeros(N,1);
errormiss = zeros(N, length(missing_ind));
% define other needed vars  
dim = size(X);
X_pred = zeros(dim(1),dim(2),N, length(conc_interval));
center = 1;
scale = 1;
modeINDAFAC =1;
method = 'Broindafac';
%ind = 0; % index for the missing data arrays
mse_threshold = 50;
% time the process 
tic 
for c = 1:length(conc_interval)
    conc = conc_interval(c);
    disp('Concentration')
    disp(conc)
    T = readtable(filename, 'Sheet', num2str(conc));
    Xs = table2array(T);
    
    filled_ind = find(~isnan(X));
    missing_ind = find(isnan(X));

    % Find the correct number of factors for this percent missing 
    parfor n = 4:N % factors (== principal components) 
        % initialise random number generator 
        disp('n')
        disp(n) 
        %perform INDAFAC with missing values on the matrix with one entry missing
        %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
        Fi = ini(X,n,1);
        [F,D, X_predloop]=missing_indafac(X,n,modeINDAFAC, center,scale,1e-10,maxiter,method, Fi);
        Xm = nmodel(F);


        % averages of metrics 
        mseloocv = (sum((X(filled_ind)-Xm(filled_ind))'.^2))/length(filled_ind);

        %if the correct starting value was not used, the error will be very
        %great 
        % can make this a for loop for future code 
        mseloop = zeros(3,1);
        randind = 1;
        mseloop(randind)=mseloocv;

        while mseloocv > mse_threshold && randind<=3
            randind=randind+1;
            rng(randind,'twister')
            %refit 
            [F,D, X_predloop]=missing_indafac(X,n,modeINDAFAC, center,scale,1e-10,maxiter,method, Fi);
            Xm = nmodel(F);

            % averages of metrics 
            mseloocv = (sum((X(filled_ind)-Xm(filled_ind))'.^2))/length(filled_ind);
            mseloop(randind)=mseloocv;
        end 
        %find random start with the lowest mse of those calculated 
        randstart = find(mseloop==min(mseloop)); 
        rng(randstart(1,1),'twister')
        [F,D, X_predloop]=missing_indafac(X,n,modeINDAFAC, center,scale,1e-10,maxiter,method, Fi);
        Xm = nmodel(F);

        errorfill(n,c, :) = (X(filled_ind)-Xm(filled_ind)); % actual prediction error
        % averages of metrics 
        msefill(n,c) = sum(errorfill(n,c, :).^2)/length(filled_ind);
       
        RAD(n,c, :) = reshape(errorfill(n,c, :), length(filled_ind),1)./(X(filled_ind));
        

        randomstart(n,c) = randstart(1,1);
        %built in metric
        %fit(n,count) = D.fit; 

        %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)

        % find the model 
        X_pred(:,:,c,n)=Xm;
        mse(n,c) = sum(errorfill(n,c,:).^2)/length(filled_ind);
        tempRAD = RAD(n,c,:);
        tempRAD = tempRAD(find(~isnan(tempRAD)));
        tempRAD = tempRAD(find(~isinf(tempRAD)));
        %RAD(n,c,1:size(tempRAD,3)) = tempRAD;
        %RAD(n,c,size(tempRAD,3):end)=0;
        RADfill(n,c)= sqrt((sum(tempRAD)).^2)/size(tempRAD,3));
        
             
    end
        
        % Find true number of factors for this % missing 
        %Find the optimal rank prediction
        [sortmse, indexsortmse(:,c)]=sort(mse(:,c)); % sorts in ascending order 

        minmse(c) = sortmse(1);
        numberOfFactors(c) = indexsortmse(1,c);
        cmin = coreconsistency(numberOfFactors(c));
        dof(c) = dim(1)*dim(2)-numberOfFactors(c)*(dim(1)+dim(2)-2);
        minaapd(c) = aapdfill(numberOfFactors(c),c);
    
end
toc % end timer

% Save variables from the run 
filenametemp = strcat('2wayrunPARAFACPoly',num2str(T),date, '.mat');
save(filenametemp)
%retrieve data using load(filename.mat)
%% Use parafac and loocv for 2-way array completion to find the error of the best rank 
%assumes no missing values in any row/column
%Find the correct number of factors 
conc_interval = 0.1:0.1:0.9;%K
%  number of factors maximum to try 
N=6;
maxiter=10000;

missing_ind = find(isnan(X));
[row,col]=find(~isnan(X));
filled_ind = find(~isnan(X));
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%

% for each concentration and possible number of factors
minmse = zeros(1,length(conc_interval));
minaapd = zeros(1,length(conc_interval)); % average absolute percent deviation = average(abs(Xpred-Xtrue)/Xtrue))
numberOfFactors = zeros(1,length(conc_interval));
dof = zeros(1,length(conc_interval));
mse = zeros(N,length(conc_interval));
msefill = zeros(N,length(conc_interval));
RADfill = zeros(N,length(conc_interval));
aapdfill = zeros(N,length(conc_interval));
indexsortmse=zeros(N,length(conc_interval));
corecon = zeros(N,length(conc_interval));
coreconsistency = zeros(N,1);
randomstart = zeros(N,length(conc_interval), length(filled_ind));

%metrics defined
smse = zeros(N,1);
fit = zeros(N,1);
it = zeros(N,1);
error = zeros(N, length(conc_interval),length(filled_ind));%missing_ind is biggest for 90% missing data, as initialisation was
errorfill = zeros(N, length(conc_interval),length(filled_ind));
RAD = zeros(N, length(conc_interval),length(filled_ind));
% metrics to compare to truth 
msemiss = zeros(N,1);
errormiss = zeros(N, length(missing_ind));
% define other needed vars  
dim = size(X);
X_pred = zeros(dim(1),dim(2),N, length(conc_interval));
center = 1;
scale = 1;
modeINDAFAC =1;
method = 'Broindafac';
%ind = 0; % index for the missing data arrays
mse_threshold = 50;

% time the process 
tic 
for c = 1:length(conc_interval)
    conc = conc_interval(c);
    disp('Concentration')
    disp(conc)
    T = readtable(filename, 'Sheet', num2str(conc));
    Xs = table2array(T);
    
    filled_ind = find(~isnan(X));
    missing_ind = find(isnan(X));

    % Find the correct number of factors for this percent missing 
    for n = N % factors (== principal components) 
        % initialise random number generator 
        disp('n')
        disp(n) 
        OldLoad = ini(X,n,2);
        [OldLoad,D, X_predloop]=missing_indafac(X,n,modeINDAFAC, center,scale,1e-3,maxiter,method, OldLoad);
        %perform INDAFAC with missing values on the matrix with one entry missing
        %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
        parfor ind = 1:length(filled_ind)
            % temporary X
            X_b = X;
            filled_index = filled_ind(ind);
            X_b(row(ind),col(ind)) = nan;% remove a true data point from Xs
            if   find(~isnan(X_b(row(ind),:)))&  find(~isnan(X_b(:,col(ind)))) & (X_b(row(ind),col(ind)))~=0 %ensure at least one value in each column and row
                 
                [F,D, X_predloop]=missing_indafac(X_b,n,modeINDAFAC, center,scale,1e-3,maxiter,method, OldLoad);
                Xm = nmodel(F);
          
                
                % averages of metrics 
                mseloocv = (sum((X(filled_ind)-Xm(filled_ind))'.^2))/length(filled_ind);

                %if the correct starting value was not used, the error will be very
                %great 
                % can make this a for loop for future code 
                mseloop = zeros(6,1);
                randind = 1;
                mseloop(randind)=mseloocv;
                
                while mseloocv > mse_threshold && randind<=6
                    randind=randind+1;
                    rng(randind,'twister')
                    %refit 
                    [F,D, X_predloop]=missing_indafac(X_b,n,modeINDAFAC, center,scale,1e-3,maxiter,method, OldLoad);
                    Xm = nmodel(F);
                   
                    % averages of metrics 
                    mseloocv = (sum((X(filled_ind)-Xm(filled_ind))'.^2))/length(filled_ind);
                    mseloop(randind)=mseloocv;
                end 
                %find random start with the lowest mse of those calculated 
                randstart = find(mseloop==min(mseloop)); 
                rng(randstart(1,1),'twister')
                [F,D, X_predloop]=missing_indafac(X_b,n,modeINDAFAC, center,scale,1e-3,maxiter,method, OldLoad);
                Xm = nmodel(F);
                
                errorfill(n,c,ind) = (X(filled_index)-Xm(filled_index)); % actual prediction error
                % averages of metrics 
                msefill(n,c,ind) = sum(errorfill(n,c,ind).^2)/length(filled_ind);
                if X(filled_index)~=0
                    RAD(n,c,ind) = errorfill(n,c,ind)/(X(filled_index));
                end 
                
                randomstart(n,c,ind) = randstart(1,1);
                %built in metric
                %fit(n,count) = D.fit; 
                error = (X(filled_ind)-Xm(filled_ind))';
                
                %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
                
                % find the model 
                X_pred(:,:,ind,n)=Xm;
            end 
        end
        mse(n,c) = sum(errorfill(n,c,:).^2)/length(filled_ind);
        tempRAD = RAD(n,c,:);
        tempRAD = tempRAD(find(~isnan(tempRAD)));
        tempRAD = tempRAD(find(~isinf(tempRAD)));
        RAD(n,c,1:size(tempRAD,3)) = tempRAD;
        RAD(n,c,size(tempRAD,3):end)=0;
        RADfill(n,c)= sqrt((sum(tempRAD.^2))/size(tempRAD,3));
    end
    % Find true number of factors for this % missing 
    %Find the optimal rank prediction
    [sortmse, indexsortmse(:,c)]=sort(mse(:,c)); % sorts in ascending order 
    
    minmse(c) = sortmse(1);
    numberOfFactors(c) = indexsortmse(1,c);
    cmin = coreconsistency(numberOfFactors(c));
    dof(c) = dim(1)*dim(2)-numberOfFactors(c)*(dim(1)+dim(2)-2);
    minaapd(c) = aapdfill(numberOfFactors(c),c);
    
end
toc % end timer

% Save variables from the run 
filenametemp = strcat('2wayrunPARAFACPoly288.15n=6',date, '.mat');
save(filenametemp)
%retrieve data using load(filename.mat)
%%
 [F,D, X_predbest]=missing_indafac(X,5,modeINDAFAC, center,scale,1e-3,1000, method);
 
 Xm=nmodel(F);
%% LOOCV for the correct number of factors 
% must be run after the correct nnumber of factors is found 
rng(2,'twister')
%Choose number of factors
numberOfFactorsLOOCV = 1;% = numberOfFactors(missing==missingLOOCV);

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
RAD = zeros(length(conc_interval),N);
mseN = zeros(1,N);
cN = zeros(1,N);
coreconsistencyN = zeros(1,N);
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%
% define how to use method again 
center = 1;
scale = 1;
modeINDAFAC =1;
method = 'Broindafac';

% define other needed vars 
X_predN = zeros(dim(1),dim(2),dim(3),N);
n = 1; % factors (== principal components) 
% LOOCV
ind=0;% counter for LOOCV
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
        ind=ind+1; % allowed, therefore the counter is increased 
        %perform INDAFAC with missing values on the matrix with one entry missing
        %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
        [F,D, X_pred(:,:,:,ind)]=missing_indafac(X_b,numberOfFactorsLOOCV,modeINDAFAC, center,scale,1e-3,1000,method);
        Xm = nmodel(F);
        %built in metric
        fitN(n,ind) = D.fit;
        %own metrics 
        % fix this 
        errorN(:,ind) = X(row(ind),col(ind),:)-Xm(row(ind),col(ind),:);
        RAD(:,ind) = errorN(:,ind)./reshape(X(row(ind),col(ind),:),size(errorN(:,ind)));
        tempvalue = 0;
        %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
        cN(n,ind)= corcond(Xm,F,[],1);
        while (cN(n,ind)<90 && errorN(n,ind)>50) && tempvalue<10
            tempvalue = tempvalue+1;
            rng(tempvalue, 'twister') % new random values 
            [F,D, X_pred(:,:,:,ind)]=missing_indafac(X_b,numberOfFactorsLOOCV,modeINDAFAC, center,scale,1e-3,1000,method);
            Xm = nmodel(F);
        end 
        
    end 
end %end LOOCV 
% averages of metrics for LOOCV 
msebest = sum(errorN(n,:).^2)./length(errorN(n,:));
RAD=RAD(find(~isinf(RAD)));
AAD_overall = sum(RAD(n,:).^2)./length(RAD(n,:));
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

function [F,D, X_pred]=missing_indafac(X,fn,mode, center,scale,conv,max_iter,method, Fi)
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
Options = zeros(6,1);
Options(6)=max_iter;
Options(2)=2;
const = 0;
OldLoad = Fi;

    dim=size(X);
    [i,j]=find(isnan(X)); % missing indices 
    if method == 'Broindafac'
        % use their indafac algorithm, not my method 
        Fi = ini(X, fn, 1); % initialise the three loadings 
        
        [F,D]=parafac(X,fn, Options, const, OldLoad);
        %function [Factors,it,err,corcondia]=parafac(X,Fac,Options,const,OldLoad,FixMode,Weights);
        Xm = nmodel(F);
        X_pred = X;
        % fill misssing values 
        for d = 1:length(i)
            X_pred(i(d),j(d))= Xm(i(d),j(d));
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

function [X_filled, missing] = filldata3(X)
    % filldata3 fills a 3-way array with the arithmetic mean of the other
    % values in the array
    % Input 
    % X = data array 
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
