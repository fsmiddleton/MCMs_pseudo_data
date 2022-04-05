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
missing =30;

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
N=6; %  number of factors 
mse = zeros(N,1);
smse = zeros(N,1);
c = zeros(N,1);
for n = 1:N
    Fi = ini(X, n, 1); % initialise the three loadings 
    [F, D]= INDAFAC(X, n, Fi); % perform INDAFAC on the matrix
    % D provides diagnostics, F is the factors 
    alphabet = 'ABCD'; %will use maximum 4way data
    for i=1:3 %each loading for each dimension of the data
        eval([alphabet(i) ,'= F{i};']);
    end 
    Xf = nmodel(F);
    error = Xtrue-Xf;
    mse(n) = sum(error.^2,'all')/(m1*m2*m3);
    smse(n) = sqrt(mse(n));
    c(n)= corcond(X,F);
end
%%
subplot(3,1,1)
plot(1:N, c,  'k.','MarkerSize',15)
ylabel('Consistency')

subplot(3,1,2)
plot(1:N, mse, 'k.','MarkerSize',15)
ylabel('Mean squared Error')

subplot(3,1,3)
plot(1:N, smse,  'k.','MarkerSize',15)
xlabel('Components')
ylabel('Error')

