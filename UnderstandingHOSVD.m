%% Trying tucker 14/11/2022
clear
clc                           
 A(:,:,1) = [.9073 .7158 -.3698; .8924 -.4894 2.4288; 2.1488 .3054 2.3753];
 A(:,:,2) = [1.7842 1.6970 .0151 ; 1.7753 -1.5077 4.0337; 4.2495 .3207 4.7146];
 A(:,:,3) = [2.1236 -0.074 1.4429; -0.6631 1.9103 -1.7495; 1.8260 2.1335 -.2716];
%A=X;
[a1, a2, a3] = nshape(A);
%%
clc
clear
clf
 load("claus.mat") %5x201x61
 X = permute(X,[2,3,1]);
% load("howto3.mat") %15x29x32
% load("howto1.mat") %10x18x5
% load("sugar.mat") %268x53x7
% load("heUNIFACforT=288.15.mat") % use multiple to form 3-way array

%load('HE4wayArrayPolySmall4', 'HE_data_sparse', 'Temps', 'conc_interval','mixture')
%A = reshape(HE_data_sparse(:,:,10 ,:),21,21,[]);
% X=A;
%fill_ind = find(~isnan(A));
%remove_ind = find(isnan(A));
perc_remove = 0.80; 
dim = size(X);
for i =1:dim(3)
    [A(:,:,i),~,~]=remove_matrix(X(:,:,i),perc_remove);
end 
miss_ind=find(isnan(A));
fill_ind = find(~isnan(A));
fns = 3;
%%
count = 0;
A(A==0)=nan;
center = 0;
scale = 0;
maxiter = 1000;
threshold = 0.6*(size(X,3)-1)/sqrt(size(X,3));
broprocess = 1;
for fn1 = fns
    count = count+1;
    fn = [fn1,fn1,4];
    Afilled = filldata3(A,'row');
    % 
    dimA=size(A);


    %A1 prediction 
    for iter =1:maxiter
        %center and scale 
        if broprocess 
            [Ac,mX,sX]=nprocess(Afilled,[1 1 1],[1 1 1]);
        else 
            [A1, ~,~]=nshape(Afilled);
            if center ==1 
                mx = mean(A1);
                A1 = A1-ones(size(A1,1),1)*mx; %columnwise mean 
            end 
            if scale ==1 
                sj = sqrt(sum((A1).^2,2));
                A1 = A1./(sj*(ones(1,size(A1,2))));
            end 
            Ac = reshape(A1,dimA);
        end 

        [A1, A2, A3]=nshape(Ac);

        %HOSVD
        %find all SVD decompositions 
        [U1,~,~]=svd(A1);
        [U2,~,~]=svd(A2);
        [U3,~,~]=svd(A3);
        %find the core tensor 
        S1 = U1'*A1*(kron(U2,U3));
        % use the number of factors you wish to use 
        U1 =  U1(:,1:fn(1));
        U2 =  U2(:,1:fn(2));
        U3 =  U3(:,1:fn(3));
        S1 = reshape(S1,dimA);
        S1pred = reshape(S1(1:fn(1),1:fn(2),1:fn(3)),fn(1),[]);
        A1_pred = U1*S1pred*(kron(U2,U3))';
        %uncentre and unscale
        if broprocess 
            %refold
            A1_pred = reshape(A1_pred,dimA);
            [A1_pred,mX,sX]=nprocess(A1_pred,[1,1,1],[1,1,1],mX,sX,-1);
        else 
            if scale ==1
                A1_pred = A1_pred.*(sj*ones(1,size(A1,2)));
            end
            if center ==1
                A1_pred = A1_pred+ones(size(A1,1),1)*mx;
            end 
            %refold
            A1_pred = reshape(A1_pred,dimA);
        end 
        Afilled(miss_ind) = A1_pred(miss_ind);
     
        if mod(iter,55)==0  %every 55 iterations, or if the convergence criterion became worse 
            [rowloop,colloop] = find(isnan(X(:,:,1))); %find missing values in this Temperature slice   
            for ind = 1:length(rowloop)
                %extract predictions for this missing value 
                tempX = Afilled(rowloop(ind),colloop(ind),:);
                if any(abs(tempX)>1200) && iter<maxiter %don't fix on the last iteration 
                    tempX(abs(tempX)>1200)=median(tempX); %set great outliers to nan to recalculate 
                    tempScore = (tempX-mean(tempX))/std(tempX); %calculate scores
                     if isnan(abs(tempScore)>threshold)
                         tempX(abs(tempScore)>threshold) = median(tempX);
                     end %remove outliers 
                    Afilled(rowloop(ind),colloop(ind),:) = tempX;
                end 
            end 
        end 
    end 
%     clear mx sj
%     
%     % A2 prediction 
%     [~,A2,~]=nshape(Afilled);
%     %Ac = reshape(A2',dimA);
%     %center and scale 
%     
%     
%     if center ==1 
%         mx = mean(A2);
%         A2 = A2-ones(size(A2,1),1)*mx; %columnwise mean 
%     end 
%     if scale ==1 
%         sj = sqrt(sum((A2).^2,2));
%         disp('scale preprocess A2')
%         disp(size(sj))
%         disp(size(A2))
%         A2 = A2./(sj*(ones(1,size(A2,2))));
%     end 
%     Ac2 = permute(reshape(A2,dimA(2),dimA(1),dimA(3)), [1 2 3]);
%     [A1, A2, A3]=nshape(Ac2);
% 
%     %HOSVD
%     %find all SVD decompositions 
%     [U1,~,~]=svd(A1);
%     [U2,~,~]=svd(A2);
%     [U3,~,~]=svd(A3);
%     %find the core tensor 
%     S2 = U2'*A2*(kron(U3,U1));
%     % use the number of factors you wish to use 
%     U1 =  U1(:,1:fn(1));
%     U2 =  U2(:,1:fn(2));
%     U3 =  U3(:,1:fn(3));
% 
%     S2 = reshape(S2,dimA(2),dimA(1),dimA(3));
%     S2pred = reshape(S2(1:fn(2),1:fn(1),1:fn(3)),fn(2),[]);
%     A2_pred = U2*S2pred*(kron(U3,U1))';
%     %uncentre and unscale
%     if scale ==1
%         disp('scale postprocess A2')
%         disp(size(A2_pred))
%         disp(size(sj))
%         A2_pred = A2_pred.*(sj*ones(1,size(A2,2)));
%     end
%     if center ==1
%         A2_pred = A2_pred+ones(size(A2,1),1)*mx;
%     end 
%     %refold
%     ans = A2_pred;
%     A2_pred = permute(reshape(A2_pred,dimA(1),dimA(3),dimA(2)),[1 3 2]);
%     clear mx sj
    %A3 prediction 
%     [~,~,A3]=nshape(Afilled);
%     %center and scale 
%     if center ==1 
%         mx = mean(A3);
%         A3 = A3-ones(size(A3,1),1)*mx; %columnwise mean 
%     end 
%     if scale ==1 
%         sj = sqrt(sum((A3).^2,2));
%         A3 = A3./(sj*(ones(1,size(A3,2))));
%     end 
%     Ac =  permute(reshape(A3,dimA(2),dimA(3),dimA(1)),[3 1 2]);
%     [A1, A2, A3]=nshape(Ac);
% 
%     %HOSVD
%     %find all SVD decompositions 
%     [U1,~,~]=svd(A1);
%     [U2,~,~]=svd(A2);
%     [U3,~,~]=svd(A3);
%     %find the core tensor 
%     S3 = U3'*A3*(kron(U1,U2));
%     % use the number of factors you wish to use 
%     U1 =  U1(:,1:fn(1));
%     U2 =  U2(:,1:fn(2));
%     U3 =  U3(:,1:fn(3));
%     S3 = reshape(S3,dimA(3),dimA(2),dimA(1));
%     S3pred = reshape(S3(1:fn(3),1:fn(2),1:fn(1)),fn(3),[]);
%     A3_pred = U3*S3pred*(kron(U1,U2))';
%     %uncentre and unscale
%     if scale ==1
%         A3_pred = A3_pred.*(sj*ones(1,size(A3,2)));
%     end
%     if center ==1
%         A3_pred = A3_pred+ones(size(A3,1),1)*mx;
%     end 
%     %refold
%     A3_pred = permute(reshape(A3_pred',dimA(2),dimA(3),dimA(1)),[3 1 2]);
% 

    error = A1_pred(fill_ind) - X(fill_ind);
    errormiss = A1_pred(miss_ind) - X(miss_ind);
    %error2 = A2_pred - Afilled;
%     error3 = A3_pred - Afilled;
    %Apred = mean(A1_pred,A2_pred,A3_pred);
    %smse2 
    C(:,:,:,1) = A1_pred;
    %C(:,:,:,2) = A2_pred;
%     C(:,:,:,2) = A3_pred;
    %M = cat(1, C{:});
    smsefill(count) = sqrt(mean(error.^2));
    wsmsefill(count) = sqrt(find_wmse_error(error(:),length(error(:))));
    smsemiss(count) = sqrt(mean(errormiss(:).^2));
    wsmsemiss(count) = sqrt(find_wmse_error(errormiss(:),length(errormiss(:))));
%     smse2(count) = sqrt(mean(mean(mean(error2.^2))));
%     smse3(count) = sqrt(mean(mean(mean(error3.^2))));
    
end

%% Plots
f1 = figure(1);
clf
plot(fns,smsefill, 'o', fns, wsmsefill, 'o')
legend('SMSE', 'wSMSE')
ylabel('Error filled entries')
xlabel('Number of factors')

f2 = figure(2);
clf
plot(fns,smsemiss, 'o', fns, wsmsemiss, 'o')
legend('SMSE', 'wSMSE')
ylabel('Error missing entries')
xlabel('Number of factors')

f3 = figure(3);
clf
plot(A(fill_ind),A1_pred(fill_ind), '.', 'MarkerSize', 9);
hold on
plot(X(miss_ind),A1_pred(miss_ind), '.', 'MarkerSize', 9);
hold on
plot([min(A(fill_ind)) max(A(fill_ind))],[min(A(fill_ind)) max(A(fill_ind))], 'k-')
hold off
legend('Filled entries', 'Missing entries')
xlabel('Original data')
ylabel('Predictions')

%% Centring each of the answers 
%this performs worse 
clear
clc
A(:,:,1) = [.9073 .7158 -.3698; .8924 -.4894 2.4288; 2.1488 .3054 2.3753];
A(:,:,2) = [1.7842 1.6970 .0151 ; 1.7753 -1.5077 4.0337; 4.2495 .3207 4.7146];
A(:,:,3) = [2.1236 -0.074 1.4429; -0.6631 1.9103 -1.7495; 1.8260 2.1335 -.2716];
%A=X;
fn = [10,10,4];
center =1;
scale =1;

load('HE4wayArrayPolySmall4', 'HE_data_sparse', 'Temps', 'conc_interval','mixture')
A = reshape(HE_data_sparse(:,:,1,:),21,21,4);
A = filldata3(A,'avr',mixture,conc_interval,'none',Temps);
dimA=size(A);
%center and scale 
A1 = reshape(A,dimA(1),[]); %mode1 unfolding 
A3 = reshape(A1',dimA(3),[]); %mode 3 unfolding 
A2 = reshape(A3', dimA(2), []); %mode 2 unfolding
if center ==1 
    A1 = A1-ones(size(A1,1),1)*mean(A1); %columnwise mean 
    A2 = A2-ones(size(A2,1),1)*mean(A2); %columnwise mean 
    A3 = A3-ones(size(A3,1),1)*mean(A3); %columnwise mean 
end 
if scale ==1 
    A1 = A1./(sqrt(sum((A1).^2,2))*(ones(1,size(A1,2))));
    A2 = A2./(sqrt(sum((A2).^2,2))*(ones(1,size(A2,2))));
    A3 = A3./(sqrt(sum((A3).^2,2))*(ones(1,size(A3,2))));
end 

%HOSVD
%find all SVD decompositions 
[U1,~,~]=svd(A1);
[U2,~,~]=svd(A2);
[U3,~,~]=svd(A3);
%find the core tensor 
S1 = U1'*A1*(kron(U2,U3));
%use the number of factors you wish to use 
U1 =  U1(:,1:fn(1));
U2 =  U2(:,1:fn(2));
U3 =  U3(:,1:fn(3));
S1 = reshape(S1,dimA);
S1pred = reshape(S1(1:fn(1),1:fn(2),1:fn(3)),fn(1),[]);
A1_pred = U1*S1pred*(kron(U2,U3))';

% post process = uncenter, unscale 
% A1 is used here for results 
if scale ==1
    A1_pred = A1_pred.*(sqrt(sum((A1).^2,2))*ones(1,size(A1,2)));
end
if center ==1
    A1_pred = A1_pred+ones(size(A1,1),1)*mean(A1);
end 

error = A1_pred - A1;
smse = sqrt(mean(mean(error.^2)))

%%
%A1 centered and scaled, reformed to A2 and A3
%A3 = reshape(A1',dimA(3),[]); %mode 3 unfolding 
%A2 = reshape(A3', dimA(2), []); %mode 2 unfolding

%%
load('HE4wayArrayPolySmall4', 'HE_data_sparse', 'Temps', 'conc_interval','mixture')
count = 0;
fns = 1:21;

for fn1 = fns
    count = count+1;
    fn = [fn1,fn1,4];
    A = reshape(HE_data_sparse(:,:,10 ,:),21,21,[]);
    A(A==0)=nan;
    center = 1;
    scale = 1;
    Afilled = filldata3(A,'reg',mixture,conc_interval,'none',Temps);
    % 
    dimA=size(A);


    %A1 prediction 
    [A1, ~,~]=nshape(Afilled);
    %center and scale 
    if center ==1 
        mx = mean(A1);
        A1 = A1-ones(size(A1,1),1)*mx; %columnwise mean 
    end 
    if scale ==1 
        sj = sqrt(sum((A1).^2,2));
        A1 = A1./(sj*(ones(1,size(A1,2))));
    end 
    Ac = reshape(A1,dimA);
    [A1, A2, A3]=nshape(Ac);

    %HOSVD
    %find all SVD decompositions 
    [U1,~,~]=svd(A1);
    [U2,~,~]=svd(A2);
    [U3,~,~]=svd(A3);
    %find the core tensor 
    S1 = U1'*A1*(kron(U2,U3));
    % use the number of factors you wish to use 
    U1 =  U1(:,1:fn(1));
    U2 =  U2(:,1:fn(2));
    U3 =  U3(:,1:fn(3));
    S1 = reshape(S1,dimA);
    S1pred = reshape(S1(1:fn(1),1:fn(2),1:fn(3)),fn(1),[]);
    A1_pred = U1*S1pred*(kron(U2,U3))';
    %uncentre and unscale
    if scale ==1
        A1_pred = A1_pred.*(sj*ones(1,size(A1,2)));
    end
    if center ==1
        A1_pred = A1_pred+ones(size(A1,1),1)*mx;
    end 
    %refold
    A1_pred = reshape(A1_pred,dimA);
    error = A1_pred - Afilled;
    %Afilled(
    smsefill(count) = sqrt(mean(mean(mean(error.^2))));
end 
for i=1:dimA(1)
    sigma1(i) = norm(reshape(S1(i,:,:),21,[]),'fro');
end 
subplot(2,1,1)
plot(fns,smsefill)
xlabel('Number of factors')
ylabel("SMSE (J/mol)")
subplot(2,1,2)
plot(1:dimA(1)-1, sigma1(1:dimA(1)-1))
xlabel('Number of factors')
ylabel('Singular value')

%% Singular values 
for i=1:dimA(1)
    sigma1(i) = norm(reshape(S1(i,:,:),21,[]),'fro');
    sigma2(i) = norm(reshape(S1(:,i,:),21,[]),'fro');
end 
for i=1:dimA(3) 
    sigma3(i) = norm(reshape(S1(:,:,i),21,[]),'fro');
end
%%
subplot(3,1,1)
plot(1:dimA(1), sigma1)
xlabel('Number of factors')
ylabel('Singular value')
subplot(3,1,2)
plot(1:21, sigma2)
xlabel('Number of factors')
ylabel('Singular value')
subplot(3,1,3)
plot(1:4, sigma3)
xlabel('Number of factors')
ylabel('Singular value')


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


