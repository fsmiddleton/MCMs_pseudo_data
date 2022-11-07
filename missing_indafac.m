function [X_pred,iter,F,err] = missing_indafac(X,fn,max_iter,conv,scale,center, fillmethod, orth,mixtures,concinterval, whichX, Temps)
   
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
    % fillmethod = method used to guess initial values of X
    % orth = number of orthogonal modes
    % mixtures = numeric respresentation of mixtures
    % concinterval = row vector of concentrations
    % whichX = indicates which form of data it is receiving, necessary for initial guesses: 'scale' or 'sign' or 'none'
    %
    % Output 
    % X_pred = filled X with new values for the missing entries
    % iter = number of iterations the algorithm used
    % F = factors used to construct the final solution
    % err = column vector of predictions
    
    dim=size(X);
    missing_ind = find(isnan(X));
    filled_ind = find(~isnan(X));
    %choose indices to use for convergence 
    indices=missing_ind;
    const=zeros(1,length(dim)); % orthogonal factors
    if orth ==1
        const(1) = 1;
    elseif orth >1
        const(1:orth)=1;
    end 
     
    if any(isnan(X)) % there is missing data
        if length(Temps)>1 % more than 3 -way data = 4-way data
            Xfilledini = zeros(dim);
            if length(dim)>3
                for i = 1:dim(4)
                    Xfilledini(:,:,:,i) = filldata3(X(:,:,:,i),fillmethod, mixtures,concinterval, whichX, Temps(i));
                end 
            else
                for i = 1:dim(3)
                    Xfilledini(:,:,i) = filldata3(X(:,:,i),fillmethod, mixtures,concinterval, whichX, Temps(i));
                end 
            end 
        else
            Xfilledini = filldata3(X, fillmethod, mixtures,concinterval, whichX, Temps); 
        end
        
        % PARAFAC options 
         Options(1) = conv;
         Options(2) = 1;
         Options(3) = 0;
         Options(4) = 0;
         Options(5) = 1000;
         Options(6) = 2000;
         
        %initialise the factors 
        DimX = size(Xfilledini);
        Xfilled = reshape(Xfilledini,DimX(1),prod(DimX(2:end)));
        
        if scale ==1 
            sj =sqrt(sum((Xfilled).^2,2)); % column vector 
            Xfilled = Xfilled./(sj*(ones(1,size(Xfilled,2))));% scale across rows 
        end 
        if center ==1  
            % recenter after scaling - optional, but should be done 
            mx2 = mean(Xfilled);
            Xfilled = Xfilled-ones(size(Xfilled,1),1)*mx2; 
        end 
        Xfilled = reshape(Xfilled,DimX);
        
        % initialises the factors 
        [Fi,~]=parafac(Xfilled,fn, Options, const);
        
        % fill predicted values into array
        %complete the array using indafac
        [F,diagnos]=INDAFAC(X,fn,Fi);
        model = nmodel(F);
        X_pred = model;
        iter = diagnos.it(2);
        % uncenter and unscale and uncenter 
        X_pred = reshape(X_pred, DimX(1), prod(DimX(2:end)));
        if center ==1
            X_pred = X_pred + ones(size(Xfilled,1),1)*mx2;
        end 
        if scale ==1
            X_pred = X_pred.*(sj*(ones(1,size(X_pred,2))));
        end
         
        X_pred = reshape(X_pred, DimX);

        err = Xfilled(filled_ind) - X_pred(filled_ind);
    else 
        % no missing data 
        disp('no missing data')
    end %end if  else 
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