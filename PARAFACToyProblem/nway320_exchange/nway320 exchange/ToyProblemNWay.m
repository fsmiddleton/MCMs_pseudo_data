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
missing =80;

if missing ==0
    filename = 'ToyProblemData3DFull.xlsx';
else 
    filename = ['ToyProblemData3D_',num2str(missing),'%missing.xlsx'];
end 

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
        T = readtable(filename, 'Sheet', num2str(i));
        X(i,:,:)=table2array(T);
    end   
end 
%% INDAFAC
N=10; %  number of factors 
mse = zeros(N,1);
smse = zeros(N,1);
c = zeros(N,1);
fit = zeros(N,1);
it = zeros(N,1);
% can also get Rho, Lambda, EV% and Max Gr
% EV must be close to 100%
for n = 1:N
    Fi = ini(X, n, 1); % initialise the three loadings 
    [F, D]= INDAFAC(X, n, Fi, diagnostics = 'off'); % perform INDAFAC on the matrix
    % D provides diagnostics, F is the factors 
    alphabet = 'ABCD'; %will use maximum 4way data
    for i=1:3 %each loading for each dimension of the data
        eval([alphabet(i) ,'= F{i};']);
    end 
    Xf = nmodel(F);
    %built in metrics
    fit(n) = D.fit;
    it(n) = D.it(2); % new computation of the Jacobian ie update accepted 
    %own metrics 
    error = Xtrue-Xf;
    mse(n) = sum(error.^2,'all')/(m1*m2*m3);
    smse(n) = sqrt(mse(n));
    c(n)= corcond(X,F);
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


%% confirm rank with plots 
% cor consistency must be close to 100%
xplot = 1:N;
subplot(3,1,1)
plot(xplot, c(xplot),  'k.','MarkerSize',15)
ylabel('Consistency')
% These will show sharp decline and stabilisation upon discovery of the
% true rank
subplot(3,1,2)
plot(xplot, mse(xplot), 'k.','MarkerSize',15)
ylabel('Mean squared Error')

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
