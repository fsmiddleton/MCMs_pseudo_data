function [X_pred,iter,F,err] = missing_parafac3(X,fn,max_iter,conv,scale,center, fillmethod, orth,mixtures,concinterval, whichX, Temp)
   
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
    elseif orth ==2
        const(1:2)=1;
    end 
     
    if any(isnan(X)) % there is missing data
        if size(dim,2)>3 % more than 3 -way data = 4-way data
            Xfilledini = zeros(dim);
            for i = 1:dim(4)
                Xfilledini(:,:,:,i) = filldata3(X(:,:,:,i),fillmethod, mixtures,concinterval, whichX, Temp(i));
            end 
        else
            Xfilledini = filldata3(X, fillmethod, mixtures,concinterval, whichX, Temp); 
        end
        f=2*conv;
        iter = 1;
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
        mx = mean(Xfilled);% you have to center, scaling is optional 
        Xfilled = Xfilled-ones(size(Xfilled,1),1)*mx;
        if scale ==1 
            sj =sqrt(sum((Xfilled).^2,2)); % column vector 
            Xfilled = Xfilled./(sj*(ones(1,size(Xfilled,2))));% scale across rows 
        end 
        if center ==1  
            % recenter after scaling - optional, but should be done 
            mx2 = mean(Xfilled);
            Xc = Xfilled-ones(size(Xfilled,1),1)*mx2; 
        end 
        Xc = reshape(Xc,DimX);
        
        % initialises the factors 
        [F,~]=parafac(Xc,fn, Options, const);
        model = nmodel(F);
        % fill predicted values into array
        X_pred = model;

        % uncenter and unscale and uncenter 
        X_pred = reshape(X_pred, DimX(1), prod(DimX(2:end)));
        if center ==1
            X_pred = X_pred + ones(size(Xfilled,1),1)*mx2;
        end 
        if scale ==1
            X_pred = X_pred.*(sj*(ones(1,size(Xfilled,2))));
        end
        X_pred = X_pred + ones(size(Xfilled,1),1)*mx; 
        X_pred = reshape(X_pred, DimX);
        Xfilledini(indices) = X_pred(indices);
        %find sum of squares 
        SS = sum(sum(sum(Xfilledini(indices).^2)));
        
        % filling algorithm repeating 
        while iter<max_iter && f>conv
            iter = iter+1;
            SSold = SS; 
            Fi = F;
            % preprocess - scale + center unfolded array  
            DimX = size(Xfilledini);
            Xfilled = reshape(Xfilledini,DimX(1),prod(DimX(2:end)));
            mx = mean(Xfilled);% you have to center, scaling is optional 
            Xfilled = Xfilled-ones(size(Xfilled,1),1)*mx;
            if scale ==1 
                sj =sqrt(sum((Xfilled).^2,2)); % column vector 
                Xfilled = Xfilled./(sj*(ones(1,size(Xfilled,2))));% scale across rows 
            end 
            if center ==1  
                mx2 = mean(Xfilled);
                Xc = Xfilled-ones(size(Xfilled,1),1)*mx2; 
            end 
            Xc = reshape(Xc,DimX);
            
            % Find the factors  
            [F,~]=parafac(Xc,fn, Options, const, Fi);
            model = nmodel(F);
            X_pred = model;
            
            %postprocess 
            % uncenter and unscale and uncenter 
            X_pred = reshape(X_pred,DimX(1),prod(DimX(2:end)));
            if center ==1
                X_pred = X_pred + ones(size(Xfilled,1),1)*mx2;
            end 
            if scale ==1
                X_pred = X_pred.*(sj*(ones(1,size(Xfilled,2))));
            end
            X_pred = X_pred + ones(size(Xfilled,1),1)*mx;
            X_pred = reshape(X_pred,DimX);
            Xfilled = X_pred;
            % fill missing values with predictions 
            Xfilledini(indices) = X_pred(indices);
            % calculate sum of squares 
            SS = sum(sum(sum(Xfilledini(indices).^2)));
            % convergence criterion 
            f = abs(SS-SSold)/(SSold);
            
        end
        err = Xfilled(filled_ind) - X_pred(filled_ind);
    else 
        % no missing data 
        disp('no missing data')
    end %end if  else 
end
