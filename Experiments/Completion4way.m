%% Toy problem for testing Nway PCA on a simulated data set 
% FS Middleton 2022/05/04
%INDAFAC code sourced from:
% Giorgio Tomasi and Rasmus Bro
%PARAFAC and missing values
%Chemometrics and Intelligent Laboratory Systems 75(2004)163-180

%%
clc
clear
 
%% Import and export data array old 
load('HEData4waySmallPoly.mat') %data is in HE_data_sparse
%mixtures is in here, as well as temps and everything else 
dim = size(HE_data_sparse);
dim1=dim(1);
dim2=dim(2);
dim3=dim(3);
dim4=dim(4);
X = HE_data_sparse;
mixtures = mixture';

%% Import data new
load('HEData4waySmallPoly.mat') %data is in HE_data_sparse
%mixtures is in here, as well as temps and everything else 
mix_original = mixture;
dim = size(HE_data_sparse);
dim1=dim(1);
dim2=dim(2);
dim3=dim(3);
dim4=dim(4);
X = HE_data_sparse(:,:,1:5,:);
conc_interval = 0.1:0.1:0.5;
mixtureT = mixture';
[comps1,~,~]=unique(mixtureT(:,[1,2]), 'rows');
[comps2,~,~]=unique(mixtureT(:,[3,4]), 'rows');
[l,Locb] = ismember(comps2,comps1,'rows');
include2 = find(Locb==0);
%Temps in the 4-way array imported here 

% all possible components in this matrix 
comps = [comps1; comps2(include2,:)];
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end 

%% LOOCV to find rank using mse  
%time the code 
tic
% declare ranks to be tested 
fns =[1:2:6,7:1:10];
concentrations=conc_interval;
Xscale = log(sign(X).*(X)).*sign(X);
Xsign = sign(X);
Xs = Xscale;
T = 298.15;
% largest length of filled_linear_ind is for 0.3
filled_ind = find(~isnan(Xs));

% vars for missing_parafac3
maxiter = 20000;
scale = 1;
center = 1;
conv = 1e-10;
winsorized_mse = 0;
fillmethod = 'avg';
orth = 1;
whichX = 'sign';

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

for c = 1:2:length(concentrations)
    conc = concentrations(c);
    %Ts = readtable(filename, 'Sheet', num2str(conc));
    %Xs = table2array(Ts);
     if strcmp(whichX, 'scale') 
        Xs = Xscale(:,:,c);
    else 
        Xs = Xsign(:,:,c);
    end 
    
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
                
                [X_pred,iters,F,err] = missing_parafac3(X_b,fn,maxiter,conv,scale,center,fillmethod,orth, mixtures,conc, whichX);
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
    filenamesave = strcat('4wayPARAFAC-All-LOOCV-X',whichX,'-maxiter=20000-T=',num2str(T),'-c=', num2str(c), '-fillmethod=',fillmethod,'-',  date, '.mat');
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

%% Find fits of the models to the data 

% declare ranks to be tested 
fns =[1:1:10, 12:2:18];
concentrations=conc_interval;
Xscale = log(sign(X).*(X)).*sign(X);
Xsign = sign(X);
%Xscale must have zeros on the diagonal too
for j = 1:dim(4)
    for i =1:dim3
        Xtemp = Xscale(:,:,i,j);
        Xtemp(1:1+size(Xtemp,1):end) = 0; % diagonals are zero
        Xscale(:,:,i,j) = Xtemp;
    end 
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
for iter = [1:2]
    %find fill method to be used
    fillmethod = fillmethods(iter);

    %intialise the metrics 
    msefill= zeros(size(fns,2),1);
    aardfill = zeros(size(fns,2),1);
    wmsefill = zeros(size(fns,2),1);

    %time the method 
    tic
    %j=0; % counter for intervals
    plotcount = 0;
     
        if strcmp(whichX, 'scale') 
            Xs = Xscale;
        else 
            Xs = Xsign;
        end 
        for j = 1:dim(4)
            Xs(:,:,:,j) = remove_nan3(Xs(:,:,:,j);
        end 
        dim = size(Xs);
        dim1 = dim(1);
        missing_ind = find(isnan(Xs));
        [row,col] = find(~isnan(Xs));
        filled_ind = find(~isnan(Xs));

            for i=1:length(fns)
                %i=i+1;
                fn = fns(i);
                disp('rank')
                disp(fn)
                %[U,D,V,X_pred]=missing_svd(X,fn,center,scale,conv,max_iter)
                %[U,D,V,St,X_pred, iters]=missing_svd(Xs,fn,1,1,1e-8,1000);
                [X_pred,iters,F,err] = missing_parafac3(Xs,fn,maxiter,conv,scale,center, fillmethod, orth,mixtures, conc, whichX);
                %X_pred is the model predictions, not only missing values are
                %filled 
                % the actual model, not just filled values - can only use errors of filled values to find model error
                Xm= X_pred;
                msefill(i) = (sum((Xm(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
                wmsefill(i) = find_wmse(Xs(filled_ind), Xm(filled_ind), length(filled_ind));
                abserrorfill = abs(X_pred(filled_ind)-Xs(filled_ind));
                Xtemp = Xs(filled_ind);
                indices_use = find(Xtemp~=0); 
                aardfill(i) = sum(abserrorfill(indices_use)./Xtemp(indices_use))/length(Xtemp(indices_use));
            end %END FNS    
    % Export results for this filling method 
    filenamesave = strcat('4wayPARAFAC-X',whichX,'-maxiter=10000-T=',num2str(T), '-fillmethod=', fillmethod ,'-orth=', num2str(orth),'-',date, '.mat');
    save(filenamesave)
end % END FILL METHODS  
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
        Options = [1e-6,1, 1, 0, 10,2500]; % convergence criterion, SVD used for initialisation, show plots,scale data, how often to show fit(iterations), max iterations]
        const = [0,0,0,0]; % no constraints 
        %[Factors,it,err,corcondia] = parafac(X,Fac,Options,const,OldLoad,FixMode,Weights);
        [F,it,err,~]=parafac(X,fn,Options,const,Fi);
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
function [F,D, X_pred]=missing_indafac(X,fn,mode, center,scale,conv,max_iter,method, mixtures,conc_interval)
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
    if strcmp(method, 'Broindafac')
        % use their indafac algorithm, not my method 
        Xf = filldata3(X,'uni',mixtures,conc_interval);
        %Fi = ini(Xf, fn, 2); % initialise the three loadings using the full matrix filled with UNIFAC guesses
        [Fi,~]=INDAFAC(Xf,fn);
        Options.convcrit.maxiter=max_iter;
        Options.diagnostics = 'off';
        Options.display = 100;
        Options.initialisation.method = 'nipals';
        Options.initialisation.maxiter = 5;
        Options.initialisation.tol = 1e-5;
        Options.convcrit.fit = 1e-10;
        Options.convcrit.grad = 1e-9;
        Options.convcrit.par = 1e-8;
        Options.convcrit.relfit = 1e-10;
        Options.weights = [];
        Options.lambdaini = 1e-1;
        Options.lambdaudpar(1)=2;
        Options.lambdaudpar(2)=1/3;
        %Fi = ini(X, fn, 1); % initialise the three loadings 
        %
        disp('I do not make it here')
        [F,D]=INDAFAC(X,fn, Fi, Options);
        Xm = nmodel(F);
        X_pred = X;
        % fill misssing values 
        for d = 1:length(i)
            X_pred(i(d),j(d),k(d))= Xm(i(d),j(d),k(d));
        end
    elseif strcmp(method, 'Broparafac')
        for i =1:dim(4)
            Xf(:,:,:,i) = filldata3(X(:,:,:,i),'avg',mixtures,conc_interval);
        end 
        %Fi = ini(Xf, fn, 2); % initialise the three loadings using the full matrix filled with UNIFAC guesses
        [Fi,~]=parafac(Xf,fn);
        Options(1) = conv;
        Options(2) = 1;
        Options(3) = 0;
        Options(4) = 0;
        Options(5) = 100;
        Options(6) = max_iter;
        const=[1 0 0]; % orthogonal factors
        [F,D] = parafac(X,fn, Options,const, Fi);
        Xm = nmodel(F);
        X_pred = X;
        % fill misssing values 
        for d = 1:length(i)
            X_pred(i(d),j(d),k(d))= Xm(i(d),j(d),k(d));
        end
    else 
        % this only uses indafac or parafac to find the loadings, with my own
        % centering and updating 
        % much longer to run algorithm 
        if length(i)>1 % there is missing data 
            Xf = filldata3(X,'uni',mixtures, conc_interval); % fill slices with averages  
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
                if strcmp(method,'Minindafac') 
                    [F,D]=INDAFAC(Xc,fn,Fi);
                    Xm = nmodel(F);
                else %minparafac
                    disp('factors')
                    disp(fn)
                    Options(1) = conv;
                    Options(2) = 1;
                    Options(3) = 0;
                    Options(4) = 0;
                    Options(5) = 10;
                    Options(6) = 100;
                    const=[0 0 0]; % no constraint on factors
                    [F,D] = parafac(Xc,fn, Options,const);
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
    % components = components ordered according to their place on the axis 
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


