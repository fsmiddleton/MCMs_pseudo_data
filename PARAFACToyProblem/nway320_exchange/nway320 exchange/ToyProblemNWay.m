%% Toy problem for testing INDAFAC on a simulated data set 
% FS Middleton 2022/03/24
%INDAFAC code sourced from:
% Giorgio Tomasi and Rasmus Bro
%PARAFAC and missing values
%Chemometrics and Intelligent Laboratory Systems 75(2004)163-180

%%
clc
clear
m1=5;
m2=30;
m3=40;
%import true data
[Xtrue, Seed, Factors] = CreaMiss(4, [m1,m2,m3], 0.01, 0, 0, 'RMV',42);

%consider missing data 
missing =50;

filename = ['ToyProblemData3D_',num2str(missing),'%missing.xlsx'];
 

export =0; %set to 0 to import 

if export ==1
    %create the matrix with missing entries 
    [X, Seed, Factors] = CreaMiss(4, [m1,m2,m3], 0.01, 0, missing/100, 'RMV',42);
    %[X,Seed,varargout] = CreaMiss(Fac,DimX,Noise,Congruence,Missing,Mode,SD)
    %export it 
    Tf = array2table(Factors);
    writetable(Tf,filename,'Sheet','Factors')

    for i = 1:m1
        % create missing matrix
        %[Xsparse,missing_ind,filled_linear_ind]=fill_matrix(X,i);
        %Tsparse = array2table(Xsparse);
        T=array2table(reshape(X(i,:,:),m2,m3));
        writetable(T, filename, 'Sheet', num2str(i))
    end 
else 
    %import 
    Factors = readtable(filename, 'Sheet', 'Factors');
    X=zeros(m1,m2,m3);
    for i =1:m1
        Ttrue = readtable('ToyProblemData3DFull.xlsx','Sheet', num2str(i));
        Xtrue(i,:,:)= table2array(Ttrue);
        T = readtable(filename, 'Sheet', num2str(i));
        X(i,:,:)=table2array(T);
    end   
end 
%% INDAFAC tensor completion 
N=10; %  number of factors maximum to try 
mse = zeros(N,1);
smse = zeros(N,1);
c = zeros(N,1);
fit = zeros(N,1);
it = zeros(N,1);
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%

%% Preprocessing 
% Fill missing data 
missing_ind = find(isnan(X));
filled_linear_ind = find(~isnan(X));
Xfilled = filldata(X);
% Unfold X into 3 different modes of arrangement (compare these) 
[x1,x2,x3]=nshape(Xfilled); 
% choose mode (2-way array) 
Xf  = x3;
% Scale X across the mode
sj = sqrt(sum(sum((Xf).^2)));
Xf = Xf/sj;
% Center X across columns 
mx = mean(Xf);
centx1 = Xf-mx*ones(size(mx,2),1);
Xc=reshape(centx1,size(X)); % reshapes back to 3-way array 

%% Tensor completion step
for n = 1:N % factors (== principal components) 
    
    % LOOCV
    k=0;% counter for LOOCV 
    for filled_ind = filled_linear_ind' %filled_linear_ind must be a row vector  
        % temporary X
        X_b = Xc;
        col = mod(filled_ind,dim2);
        if col ==0
            col=dim2;% mod(integer*m,m)=0, therefore the column should be column m
        end 
        row = ceil(filled_ind/dim2);
        X_b(filled_ind) = nan;% remove a point from Xs
        if find(X_b(:,col)) &&  find(X_b(row,:)) %ensure at least one value in each column and row
            k=k+1;% allowed, therefore the counter is increased 
            boot_removed_col(k) = col;
            boot_removed_row(k) = row;
            %perform iterative PCA on the matrix with one entry missing
            [U,D,V,X_pred_boot(:,:)]=missing_svd(X_b,min_fn,1,1,1e-3,1000);
            %save error for this value left out 
            error(k) = X_pred_boot(filled_ind)- Xs(filled_ind);
        end 
    end %end LOOCV 
                
    Fi = ini(Xc, n, 1); % initialise the three loadings 
    [F, D]= INDAFAC(Xc, n, Fi, diagnostics = 'off'); % perform INDAFAC on the matrix
    % D provides diagnostics, F is the factors 
    alphabet = 'ABCD'; %will use maximum 4way data
    for i=1:3 %each loading for each dimension of the data
        eval([alphabet(i) ,'= F{i};']);
    end 
    Xm = nmodel(F);
    %built in metrics
    fit(n) = D.fit;
    it(n) = D.it(2); % new computation of the Jacobian ie update accepted 
    %own metrics 
    error = Xtrue-Xm;
    mse(n) = sum(error.^2,'all')/(m1*m2*m3);
    smse(n) = sqrt(mse(n));
    %[Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot)
    c(n)= corcond(Xm,F,[],1);
end
%Find the optimal rank prediction
minmse = min(mse);
% assumes that minimum rank tested was 1 
numberOfFactors = find(mse==minmse);
minsmse = smse(numberOfFactors);
cmin = c(numberOfFactors);
dof = m1*m2*m3-numberOfFactors*(m1+m2+m3-2);
%create model
Fi = ini(X, numberOfFactors, 1);
[F, D]= INDAFAC(X, numberOfFactors, Fi);
for i=1:3 %each loading for each dimension of the data
        eval([alphabet(i) ,'= F{i};']);
end 
Xm = nmodel(F);
%% Post processing 
X_pred = Xm;

[x1,x2,x3]=nshape(Xm); % unfolds X into 3 different modes of arrangement (compare these) 
% choose mode (2-way array) 
Xg  = x3;
%Un-center
Xg = Xg+mx*ones(size(mx,2),1);
% Un-scale
Xg =Xg*sj;
X_pred=reshape(Xg,size(X)); % reshapes back to 3-way array


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
function [X_filled] = filldata(X)
    % fill in slices 
    [i, j, k]=find(isnan(X));% returns rows and columns with nonzero elements
    missing_ind = find(isnan(X));% column vector 
    X_filled=X;
    X_filled(missing_ind)=0; %fill NaN values with 0
    dim = size(X);
    for d = 1:dim(3)
        Xtemp = X(:,:,d);
        mean = sum(Xtemp)./(ones(1,dim(2))*dim(3)-sum(missing_ind,1)); % (j and k)
        missing = find(isnan(Xtemp));
        [i,j]=find(isnan(Xtemp));
        mean_col = sum(Xtemp,1)./(ones(1,dim(2))*dim(1)-sum(missing,1)); %columns are dimension 1
        mean_row= sum(Xtemp,2)./(ones(1,dim(1))*dim(2)-sum(missing,2)); % rows are dimension 2 
        for k =1:length(i) 
           % for all NaN elements that exist, loop through them to replace
            % with means of the slice
            X_filled(i(k),j(k), d)=(mean_row(i(k))+mean_col(j(k)))/2;
        end 
    end     
end 
