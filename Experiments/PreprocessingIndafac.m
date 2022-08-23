
%% Preprocessing development 
% Option 1 = copy PCA template which uses missing_svd 
% INDAFAC is applied to the filled matrix iteratively and is only used to find the factors. 
%% Import toy data 
clc 
clear

missing =50;
filename = ['ToyProblemData3D_',num2str(missing),'%missing.xlsx'];

dim1 = 5;
dim2=30;
dim3=40;
%import 

Factors = readtable(filename, 'Sheet', 'Factors');
X=zeros(dim1,dim2,dim3);
for i =1:dim1
    Ttrue = readtable('ToyProblemData3DFull.xlsx','Sheet', num2str(i));
    Xtrue(i,:,:)= table2array(Ttrue);
    T = readtable(filename, 'Sheet', num2str(i));
    X(i,:,:)=table2array(T);
end 

%% Preprocessing steps 
[Xf, miss]=filldata3(X);

%%
[i, j, k]=findnan3(X);% returns rows and columns with nonzero elements
    dim = size(X);
    missing = isnan(X); % 
    X_filled = X;
    X_filled(missing)=0;
%%
    % fill missing values with zeros 
%     for ind =1:length(i)
%         X_filled(i(ind),j(ind),k(ind))=0;
%     end
    % find all means 
    mean1 = sum(X_filled,1)./(ones(1,dim(2),dim(3))*dim(1)-sum(missing,1)); % (1,j,k) dim1 
    mean2= sum(X_filled,2)./(ones(dim(1),1,dim(3))*dim(2)-sum(missing,2)); % (i,1,k)  dim2 
    mean3 = sum(X_filled,3)./(ones(dim(1),dim(2),1)*dim(3)-sum(missing,3));% (i,j,1) dim3 
    %replace mean = nan with mean =0 
    mean1(find(isnan(mean1)))=0;
    
    disp(find(isnan(mean1)))
    disp(find(isnan(mean2)))
    disp(find(isnan(mean3)))
    for ind =1:length(i)
       % for all NaN elements that exist, loop through them to replace
        X_filled(i(ind),j(ind), k(ind))=(mean1(1,j(ind), k(ind))+mean2(i(ind),1, k(ind))+mean3(i(ind),j(ind),1))/3;
    end      

%%

%% Functions
function [X_filled, missing] = filldata3(X)
% Input 
% X = data array 
% Output 
% X_filled = filled data array
% missing = M matrix with ones for nan values 

% filldata3 fills a 3-way array with the arithmetic mean of the other
% values in the array
    % fill as done for PCA with SVD (averages of each dimension)
    [i, j, k]=findnan3(X);% returns rows and columns with nonzero elements
    dim = size(X);
    missing = isnan(X); % 
    X_filled = X;
    X_filled(missing)=0;
    % fill missing values with zeros 
%     for ind =1:length(i)
%         X_filled(i(ind),j(ind),k(ind))=0;
%     end
    % find all means 
    mean1 = sum(X_filled,1)./(ones(1,dim(2),dim(3))*dim(1)-sum(missing,1)); % (1,j,k) dim1 
    mean2= sum(X_filled,2)./(ones(dim(1),1,dim(3))*dim(2)-sum(missing,2)); % (i,1,k)  dim2 
    mean3 = sum(X_filled,3)./(ones(dim(1),dim(2),1)*dim(3)-sum(missing,3));% (i,j,1) dim3 
    mean1(find(isnan(mean1)))=0;
    mean2(find(isnan(mean2)))=0;
    mean3(find(isnan(mean3)))=0;

    for ind =1:length(i)
       % for all NaN elements that exist, loop through them to replace
        X_filled(i(ind),j(ind), k(ind))=(mean1(1,j(ind), k(ind))+mean2(i(ind),1, k(ind))+mean3(i(ind),j(ind),1))/3;
    end      
end 



function [i,j,k]=findnan3(X)
% Input 
% X = data array 
% Output 
% i, j, k = indexes of nan values in X

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
        [itemp, jtemp]= find(isnan(Xtemp));
        i = [i; itemp];
        j = [j;jtemp];
        ktemp = ones(length(itemp),1)*d;
        k = [k;ktemp];
    end 
end

function [F,D, X_pred]=missing_indafac(X,fn,mode, center,scale,conv,max_iter)
% Input 
% X = data array 
% fn = number of factors 
% mode = mode of unfolding used 
% center  =1 center data 
% scale  =1 scale data 
%         =0 do not scale data  
% conv = convergence criterion = difference between new and old 
 
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
    [i,j,k] = findnan3(isnan(X));
    if any(isnan(X)) % there is missing data 
        Xf = filldata3(X); % fill slices with averages  
        SS = sumsquare3(Xf,i,j,k);
        
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
                    Xf = Xp/sj;
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
            % INDAFAC step 
            Fi = ini(Xc, fn, 1); % initialise the three loadings 
            [F,D]=indafac(Xc,fn,Fi);
            Xm = nmodel(F);
            
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
            SS = sumsquare3(Xf,i,j,k);
            f = abs(SS-SSold)/(SSold);
            iter = iter+1;
        end 
    else 
        % no missing data
        disp('No missing data')
    end %end if  else 
end % end function 





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