function [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd(X,fn,center,scale,conv,max_iter, use_missing,fillmethod,mixtures,whichX,conc,T)
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
    % use_missing = use the missing entries for the convergence, =1 to use
    % missing 
    % Output 
    % S,V,D from X = SVD'
    % St = SV
    % X_pred = filled X with new values for the missing entries
    [m,n]=size(X);
    missing_ind = find(isnan(X));
    filled_ind = find(~isnan(X));
    %choose indices to use for convergence 
    if use_missing ==1
        indices=missing_ind;
    else 
        indices = filled_ind;
    end 
     
    if any(isnan(X)) % there is missing data 
        %Xfilled = fill_data(X); 
        Xfilled = filldata3(X,fillmethod,mixtures,conc, whichX, T);
        SS = sum(sum(Xfilled(indices).^2));
        f=2*conv;
        iter = 0;
        while iter<max_iter && f>conv
            iter = iter+1;
            SSold = SS; 
            % preprocess = scale and then center 
            mx = mean(Xfilled);
            if scale ==1 
                sj =sqrt(sum((Xfilled).^2,2));
                Xfilled = Xfilled./(sj*(ones(1,size(Xfilled,2))));
            end 
            if center ==1
                mx = mean(Xfilled); %columnwise mean 
                Xc = Xfilled-ones(m,1)*mx; %centering of the data done - subtract mean from each entry in a column 
            else 
                Xc=Xfilled;
            end 
            %perform SVD
            [S,V,D]=svd(Xc);
            AllSVs = V;
            St=S*V;
             St = St(:,1:fn);
             V=V(:,1:fn);
             S=S(:,1:fn);
             D=D(:,1:fn);
            % post process = uncenter, unscale  
            if center ==1
                X_pred = St*D'+ones(m,1)*mx;
            else 
                X_pred = St*D';
            end 
            if scale ==1
                X_pred = X_pred.*(sj*(ones(1,size(Xfilled,2))));
            end
            
            % fill misssing values 
            Xfilled(missing_ind)=X_pred(missing_ind);
            SS = sum(sum((X_pred(indices)).^2));
   
            f = abs(SS-SSold)/(SSold);
            disp('Iter')
            disp(iter)
            disp('Loss function')
            disp(f)
        end
        iterations = iter;
        
    else 
        iterations =1;
        % no missing data 
        Xfilled=X;
        if scale ==1 
            sj = sqrt(sum(sum((Xfilled).^2)));
            Xfilled = Xfilled/sj;
        end
        if center ==1
            mx = mean(Xfilled);
            %Xc = normalize(Xf); %normalizing 
            Xc = Xfilled-ones(m,1)*mx; %centering
        else 
            Xc=X;
        end 
        [S,V,D]=svd(Xc);
        St=S*V;
        AllSVs=V;
        St = St(:,1:fn);
        V=V(:,1:fn);
        S=S(:,1:fn);
        D=D(:,1:fn);
         
        if center ==1
            X_pred = St*D'+ones(m,1)*mx;
            
        else 
            X_pred = St*D';
        end 
        if scale ==1
                X_pred = X_pred*sj;
        end
    end %end if  else 
end % end fu