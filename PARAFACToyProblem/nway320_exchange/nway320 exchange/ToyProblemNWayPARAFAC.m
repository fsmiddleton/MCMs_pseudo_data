%% Toy problem for testing Nway PCA on a simulated data set 
% FS Middleton 2022/05/04
%INDAFAC code sourced from:
% Giorgio Tomasi and Rasmus Bro
%PARAFAC and missing values
%Chemometrics and Intelligent Laboratory Systems 75(2004)163-180

%%
clc
clear
%% Import and export toy data arrays 
%dimensions of the toy array 
dim1 = 30;
dim2 = 40;
dim3 = 5;
dim4 = 7;
dim = [dim1, dim2, dim3, dim4];
%[X,SeedOut,Factors{1:length(dimX)}] = CreaMiss(Rank,dimX,Noise,Congruence,Missing,Mode,SeedIn);
%create true data 
[Xtrue, Seed, Factors] = CreaMiss(5, [dim1,dim2,dim3,dim4], 0.01, 0, 0, 'RMV',42);


missing =20;
filename = ['ToyProblemData4D_',num2str(missing),'%missing_2.xlsx'];

export =0; %set to 0 to import 

if export ==1
    %create the matrix with missing entries 
    [X, Seed, Factors] = CreaMiss(5, [dim1,dim2,dim3, dim4], 0.01, 0, missing/100, 'RMV',42);
    %[X,Seed,varargout] = CreaMiss(Fac,DimX,Noise,Congruence,Missing,modeINDAFAC,SD)
    %export it 
    Tf = array2table(Factors);
    writetable(Tf,filename,'Sheet','Factors')

    for i = 1:dim3
        for j =1:dim4
        % create missing matrix
        T=array2table(reshape(X(:,:,i,j),dim1,dim2));
        sheetname = strcat(num2str(i),';',num2str(j));
        writetable(T, filename, 'Sheet', sheetname)
        end
    end 
else 
    %import 
    missing = 50;
    filename = ['ToyProblemData4D_',num2str(missing),'%missing_2.xlsx'];
    truefilename = ['ToyProblemData4D_0%missing_2.xlsx'];
    [X, Xtrue, Factors] = importX(filename, dim, truefilename);  
end 

%% Tensor completion loops
% Tensor completion, allowing the global minimum for a number of factors
% and an array with a %missing from the true array 
clc
clear
% import data to initialise vars 
dim = [30,40,5,7];
missing = 80; % maximum amount missing -> smallest number of entries in filled_linear_ind and greatest of missing_ind
filename = ['ToyProblemData4D_',num2str(missing),'%missing_2.xlsx'];
truefilename = ['ToyProblemData4D_0%missing_2.xlsx'];
[X, Xtrue, Factors] = importX(filename, dim, truefilename);
missing_ind = find(isnan(X));
filled_linear_ind = find(~isnan(X));
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%

%Find the correct number of factors 
missing = 30:10:80;
%  number of factors maximum to try 
N=7;

% for each % missing 
minmse = zeros(1,length(missing));
numberOfFactors = zeros(1,length(missing));
dof = zeros(1,length(missing));
mse = zeros(N,length(missing));
msefill = zeros(N,length(missing));
indexsortmse=zeros(N,length(missing));
c = zeros(N,length(missing));
coreconsistency = zeros(N,1);
randomstart = zeros(N,length(missing));
iterations = zeros(N,length(missing));

%metrics defined
smse = zeros(N,1);
fit = zeros(N,1);
it = zeros(N,1);
error = zeros(N, length(missing_ind));%missing_ind is biggest for 90% missing data, as initialisation was
errorfill = zeros(N,length(missing_ind));
% metrics to compare to truth 
msemiss = zeros(N,1);
errormiss = zeros(N, length(missing_ind));
% define other needed vars  
X_pred = zeros(dim(1),dim(2),dim(3),dim(4),length(missing));
center = 1;
scale = 1;
model =1;
method = 'Broparafac';
alphabet = 'ABCD'; %will use maximum 4way data
count = 0; % index for the missing data arrays
mse_threshold = 25;
% time the code 
tic 
for miss = missing
    disp('% missing')
    disp(miss)
    count = count+1;
    %import wanted data 
    [X, Xtrue, Factors] = importX(filename, dim, truefilename);
    filled_linear_ind = find(~isnan(X));
    missing_ind = find(isnan(X));
    no_filled = miss/100*dim(1)*dim(2)*dim(3)*dim(4);
    %remove nan slabs from the 4-way array 
    [Xnew,dim, nan_coords]=remove_nan4(X);
    
    % Find the correct number of factors for this percent missing 
    for n = 1:N % factors (== principal components) 
        % initialise random number generator 
        randind = 2;
        rng(randind, 'twister')
        
        disp('n')
        disp(n) 
        %perform INDAFAC with missing values on the matrix with one entry missing
        %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
        [F,D, X_pred]=missing_parafac(X,n,model, center,scale,1e-3,1000,method);
        Xm = nmodel(F);
        errorfill(n,1:length(filled_linear_ind)) = (X(filled_linear_ind)-Xm(filled_linear_ind))'; % actual prediction error
        % averages of metrics 
        msefill(n,count) = sum(errorfill(n,1:length(filled_linear_ind)).^2)/length(filled_linear_ind);

        %if the correct starting value was not used, the error will be very
        %great 
        % can make this a for loop for future code 
        while msefill(n,count) > mse_threshold
            randind=randind+1;
            disp('randind')
            disp(randind)
            rng(randind,'twister')
            %refit 
            [F,D, X_pred]=missing_parafac(X,n,model, center,scale,1e-3,1000,method);
            Xm = nmodel(F);
            errorfill(n,1:length(filled_linear_ind)) = (X(filled_linear_ind)-Xm(filled_linear_ind))'; % actual prediction error 
            msefill(n,count) = sum(errorfill(n,1:length(filled_linear_ind)).^2)/length(filled_linear_ind);
        end 
        iterations(n,count) = D(1);
        randomstart(n,count) = randind;
        %built in metric
        error(n,1:length(filled_linear_ind)) = (Xtrue(filled_linear_ind)-Xm(filled_linear_ind))';
        mse(n,count) = sum(error(n,1:length(filled_linear_ind)).^2)/length(filled_linear_ind);
        errormiss(n,1:length(missing_ind))= (Xtrue(missing_ind)-Xm(missing_ind))';
        msemiss(n,count) = sum(errormiss(n,1:length(missing_ind)).^2)/length(missing_ind);
        %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
        c(n,count)= corcond(Xm,F,[],1);
    end
    % Find true number of factors for this % missing 
    %Find the optimal rank prediction
    [sortmse, indexsortmse(:,count)]=sort(mse(:,count)); % sorts in ascending order 
    m = 0; % counter for finding consisent number of factors with lowest mse 
    check =1;
    coreconsistency = c(:,count);
    while check==1
        m=m+1;
        if coreconsistency(indexsortmse(m,count))>0.9 % should be consistent 
            check =0;
        end 
    end 
    minmse(count) = sortmse(m);
    numberOfFactors(count) = indexsortmse(m,count);
    cmin(count) = coreconsistency(numberOfFactors(count));
    dof(count) = dim(1)*dim(2)*dim(3)-numberOfFactors(count)*(dim(1)+dim(2)+dim(3)-2);
    msemiss_min(count) = msemiss(numberOfFactors(count),count);
    % find the model 
    [F,D, X_pred(:,:,:,:,count)]=missing_parafac(X,numberOfFactors(count),model, center,scale,1e-3,1000, method);
end
% Create table with results 
Results = table(["% Missing"; (missing')],[ "Factors";numberOfFactors'],[ "MSE";minmse'],[ "Core consistency";cmin'], ["Degrees of freedom";dof'], ["MSE missing data";msemiss_min']);
disp(Results)
toc % end timer 


%% LOOCV for the correct number of factors for a % missing 
% Fix now for 4D
% must be run after the correct nnumber of factors is found 
missingLOOCV = 50;
dim = [30,40,5,7];
rng(2,'twister')
numberOfFactorsLOOCV =4;% = numberOfFactors(missing==missingLOOCV);
filename = ['ToyProblemData4D_',num2str(missingLOOCV),'%missing.xlsx'];
truefilename = ['ToyProblemData4D_0%missing.xlsx'];
[X, Xtrue, Factors] = importX(filename, dim, truefilename );
filled_linear_ind = find(~isnan(X));
missing_ind = find(isnan(X));
% find filled indices 
[i,j,k]=findfill4(X);
% define metrics here that need to be used 
N = length(filled_linear_ind);
smseN = zeros(1,N);
fitN = zeros(1,N);
itN = zeros(1,N);
errorN = zeros(1,N);
mseN = zeros(1,N);
cN = zeros(1,N);
coreconsistencyN = zeros(1,N);
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%
% define how to use method again 
center = 1;
scale = 1;
model =1;
method = 'Broparafac';

% define other needed vars 
X_predN = zeros(dim(1),dim(2),dim(3),N);
n = 1; % factors (== principal components) 
% LOOCV
count=0;% counter for LOOCV
ind = 0; % counter for filled indices 

% time the LOOCV
tic
for filled_ind = filled_linear_ind' %filled_linear_ind must be a row vector  
    % temporary X
    X_b = X;
    ind=ind+1;
    disp(ind)
    X_b(filled_ind) = nan;% remove a point from Xs
    if find(~isnan(X_b(i(ind),j(ind),:))) &  find(~isnan(X_b(i(ind),:,k(ind))))&  find(~isnan(X_b(i(ind),:,k(ind))))%ensure at least one value in each column and row
        count=count+1; % allowed, therefore the counter is increased 
        %perform INDAFAC with missing values on the matrix with one entry missing
        %[F,D, X_pred]=missing_indafac(X,fn,modeINDAFAC, center,scale,conv,max_iter)
        [F,D, X_pred(:,:,:,count)]=missing_parafac(X_b,numberOfFactorsLOOCV,model, center,scale,1e-3,1000,method);
        Xm = nmodel(F);
        %built in metric
        fit(n,count) = D.fit;
        %own metrics 
        % fix this 
        errorN(n,count) = Xtrue(filled_ind)-Xm(filled_ind);
        while errorN(n,count) >1
            rand(erorrN(n,count),'twister') % new random values 
            [F,D, X_pred(:,:,:,count)]=missing_parafac(X_b,numberOfFactorsLOOCV,model, center,scale,1e-3,1000,method);
            Xm = nmodel(F);
        end 
        %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
        cN(n,count)= corcond(Xm,F,[],1);
    end 
end %end LOOCV 
% averages of metrics for LOOCV 
msebest = sum(errorN(n,:).^2)./length(errorN(n,:));
coreconsistencybest = sqrt(sum(cN(n,:).^2)./length(cN(n,:)));

toc 
%% confirm rank with plots 
% cor consistency must be close to 100%
xplot = 1:N;
subplot(3,1,1)
plot(xplot, c(xplot),  'k.','MarkerSize',15)
ylabel('Consistency')
% These will show sharp decline and stabilisation upon discovery of the
% true rank
subplot(3,1,2)
plot(xplot, smse(xplot), 'k.','MarkerSize',15)
ylabel('Square root of the mean squared Error')

subplot(3,1,3)
plot(xplot, fit(xplot),  'k.','MarkerSize',15)
xlabel('Components')
ylabel('Fit: Loss function ')
%% Plt the missing data structure 
% these must have the same sizes as x
v=X;

xslice = [15,25];    % location of y-z planes
yslice = [2,3,4,5];              % location of x-z plane
zslice = [1,25];         % location of x-y planes
clf
slice(v,xslice,yslice,zslice)
xlabel('x')
ylabel('y')
zlabel('z')
%% 
function [X,dim, nancoords]=remove_nan4(X)
    % saves the columns and rows with only Nan values and removes them from
    % the 4-way array 
    % Input 
    % X = 4-way array 
    % Output
    % X = 4-way array without slabs of only nan values 
    % dim = dimensions of X outputted 
    % nancoords = indices of each dimension that contained only nan values.
    % Stored as a dictionary

    dim = size(X);
    % z slabs
    t1 = all(isnan(X),[1,2,3]);% finds slabs with only nan values 
    t1 = reshape(t1,[1,dim(4)]); %logical 
    r1=find(t1);
    X(:,:,:,r1)=[];

    % x slabs 
    t2 = all(isnan(X),[2,3,4]);
    t2 = reshape(t2,[1,dim(1)]); %logical 
    r2=find(t2);
    X(r2,:,:,:)=[];

    % y slabs 
    t3 = all(isnan(X),[1,3,4]);
    t3 = reshape(t3,[1,dim(2)]); %logical 
    r3=find(t3);
    X(:,r3,:,:)=[];
    
    % last slabs 
    t4 = all(isnan(X),[1,2,4]);
    t4 = reshape(t4,[1,dim(3)]); %logical 
    r4=find(t4);
    X(:,:,r4,:)=[];

     %new X 
    dim = size(X);
    nancoords{1} = r1; 
    nancoords{2} = r2; 
    nancoords{3} = r3; 
    nancoords{4} = r4; 

end 

function [X, Xtrue, Factors] = importX(filename, dim, truefilename)
% Imports the toy data stored in a spreadsheet with each matrix on its own
% sheet
% Input 
% filename = name of the file containing desired data 
% dim = dimensions of the data expected in the spreadsheet or amount to
% extract 
% truefilename = name of the file containing desired data in its
% non-noisy form without missing data 
% Output 
% X = desired 4-way array 
% Xtrue = true, non-noisy, filled, 4-way array 
% Factors = Factors used to construct the array; Xtrue = ABCD

    Factors = readtable(filename, 'Sheet', 'Factors');
    X=zeros(dim(1),dim(2),dim(3),dim(4));
    Xtrue=zeros(dim(1),dim(2),dim(3),dim(4));
    for i =1:dim(3)
        for j = 1:dim(4)
            sheetname = strcat(num2str(i),';', num2str(j));
            Ttrue = readtable(truefilename,'Sheet', sheetname);
            Xtrue(:,:,i,j)= table2array(Ttrue);
            T = readtable(filename, 'Sheet', sheetname);
            X(:,:,i,j)=table2array(T);
        end 
    end
end

function [F,Diagnostics, X_pred]=missing_parafac(X,fn,mode, center,scale,conv,max_iter,method)
% Input 
% X = data array 
% fn = number of factors 
% mode = mode of unfolding used 
% center  = 1 center data 
% scale  = 1 scale data 
%         =0 do not scale data  
% conv = convergence criterion = difference between new and old 
% method ='parafac' (own algorithm) or  'Broparafac'
 
% Output 
% F = factors (A,B,C,D)
% Diagnostics = diagnostics 
%     D.it = iterations 
%     D.fit = fit
%     D.error = error
%     D.cor = core consistency 
%     D.stop = reason for stopping (max iterations or conv reached)
% X_pred = predicted array 

% PARAFAC is used to find the factors of a 4-way array for a given number
% of factors (fn), with options as to how to unfold the matrix for centering and scaling, 
% as well as whether to center and scale and when to stop iterations 


    dim=size(X);
    missing_ind =find(isnan(X)); % missing indices 
    if strcmp(method, 'Broparafac')
        % use their parafac algorithm, not my method, which handles missing
        % values as nan values 
        Fi = ini(X, fn, 1); % initialise the three loadings 
        %[Factors,it,err,corcondia] = parafac(X,Fac,Options,const,OldLoad,FixMode,Weights);
        [F,it,err,~]=parafac(X,fn);
        Diagnostics = [it,err];
        Xm = nmodel(F);
        X_pred = X;
        % fill misssing values 
        for d = 1:length(missing_ind)
            X_pred(d)= Xm(d);
        end
    else 
        % this only uses indafac to find the loadings, with my own
        % centering and updating 
        % much longer to run algorithm 
        if length(i)>1 % there is missing data 
            Xf = filldata4(X); % fill slices with averages  
            SS = sumsquare3(Xf,i,j,k); % sum of squares of missing values 
            X_pred = X;
            f=2*conv;
            iter = 1;
            while iter<max_iter && f>conv
                SSold = SS;
                % preprocess = scale and then center 

                if scale ==1 || center ==1
                    [x1,x2,x3,x4]=nshape(Xf);
                    %For an I x J x K x L four-way array this means X1 is I x JKL, X2 is
                    % J x IKL, X3 is K x IJL, and X4 is L x IJK
                    %unfold to matrix
                    if mode ==1
                        Xp=x1;
                    elseif mode ==2
                        Xp = x2;
                    elseif mode ==3 
                        Xp = x3;
                    else
                        Xp=x4;
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

                for alph=1:4 %each loading for each dimension of the data
                    eval([alphabet(alph) ,'= Fi{alph};']);
                end 
                % indafac or parafac step 
                if strcmp(method,'parafac') 

                    [F,Diagnostics]=parafac(Xc,fn,Fi);
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
                for d = 1:length(missing_ind)
                    X_pred(d)= Xm(d);
                end 
                % check for convergence 
                % fix this function 
                SS = sumsquare4(Xf,missing); % sum of squares of missing values 
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

function [X_filled, missing] = filldata4(X)
% Input 
% X = data array 
% Output 
% X_filled = filled data array
% missing = M matrix with ones for nan values 

% filldata3 fills a 3-way array with the arithmetic mean of the other
% values in the array
    % fill as done for PCA with SVD (averages of each dimension)
    [i, j, k, m]=findnan4(X);% returns rows and columns with nan elements
    dim = size(X);
    missing = isnan(X); %linear indices 

    X_filled = X;
    % fill missing values with zeros 
    for ind = 1:length(i)
        X_filled(i(ind),j(ind),k(ind),m(ind))=0;
    end 

    % find all means 
    mean1 = sum(X_filled,1)./(ones(1,dim(2),dim(3),dim(4))*dim(1)-sum(missing,1)); % (1,j,k) dim1 
    mean2= sum(X_filled,2)./(ones(dim(1),1,dim(3), dim(4))*dim(2)-sum(missing,2)); % (i,1,k)  dim2 
    mean3 = sum(X_filled,3)./(ones(dim(1),dim(2),1, dim(4))*dim(3)-sum(missing,3));% (i,j,1) dim3 
    mean4 = sum(X_filled,4)./(ones(dim(1),dim(2),dim(3),1)*dim(4)-sum(missing,4));% (i,j,1) dim3
    %replace nan means with 0 
    mean1(find(isnan(mean1)))=0;
    mean2(find(isnan(mean2)))=0;
    mean3(find(isnan(mean3)))=0; 
    mean4(find(isnan(mean4)))=0;
    for ind =1:length(i)
       % for all NaN elements that exist, loop through them to replace
        X_filled(i(ind),j(ind), k(ind), m(ind))=(mean1(1,j(ind), k(ind),m(ind))+mean2(i(ind),1, k(ind), m(ind))+mean3(i(ind),j(ind),1, m(ind))+mean4(i(ind),j(ind),k(ind), 1))/4;
    end      
end 
function [i,j,k,m]=findfill4(X)
% Input 
% X = data array 
% Output 
% i, j, k, m = indexes of existing values in X

% findnan4 finds the nan indices of a 4-way array 
    dim3 = size(X,3);
    dim4 = size(X,4);
    i=[];
    j=[];
    k=[];
    m = [];
    for d = 1:dim3 
        for f = 1:dim4
            Xtemp = X(:,:,d,f);
            [itemp, jtemp]= find(~isnan(reshape(Xtemp,[dim1,dim2])));
            i = [i; itemp];
            j = [j; jtemp];
            ktemp = ones(length(itemp),1)*d;
            mtemp = ones(length(itemp),1)*f;
            k = [k;ktemp];
            m = [m;mtemp];
        end 
    end 
end

function [i,j,k,m]=findnan4(X)
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
    dim4 = dim(4);

    i=[];
    j=[];
    k=[];
    m=[];
    for d = 1:dim3
        for f = 1:dim4
        % done per z slice 
            Xtemp = X(:,:,d,f);
            [itemp, jtemp]= find(isnan(reshape(Xtemp,[dim1,dim2])));
            i = [i; itemp];
            j = [j; jtemp];
            ktemp = ones(length(itemp),1)*d;
            mtemp = ones(length(itemp),1)*f;
            k = [k;ktemp];
            m = [m;mtemp];
        end 
    end 
end

function [SS]=sumsquare4(X,indices)
% Input 
% X = data array 
% indices = indices of values of interest
% Output 
% SS = sum of the squares of the values of interest 
    tempX = X(indices); % linear 
    SS = sum(tempX.^2);
end 
