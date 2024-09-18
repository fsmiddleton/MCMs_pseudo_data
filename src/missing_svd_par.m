function [S,V,D,St,X_predall, AllSVs]=missing_svd_par(X,fn,center,scale,conv,max_iter, use_missing,fillmethod,mixtures,whichX,conc,T,thresholdperc)
    % Fill a matrix of missing data using PCA with SVD and a given number of
    % PCs. Can also handle non-missing data. Missing data is handled as NaN
    % values 
    %
    % Input 
    % X = 3-way array 
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
    [m,n]=size(X(:,:,1));
    missing_ind = find(isnan(X(:,:,1)));
    [row,col] =  find(isnan(X(:,:,1)));
    filled_ind = find(~isnan(X(:,:,1)));
    threshold = thresholdperc/100*(size(X,3)-1)/sqrt(size(X,3));
    %choose indices to use for convergence 
    if use_missing ==1
        indices=missing_ind;
    else %uses the filled indices to provide predictions 
        indices = filled_ind;
    end 
    SS = zeros(length(conc),1); 
    
    if any(isnan(X(:,:,1))) % there is missing data 
        Xfilledall = filldata3(X,fillmethod,mixtures,conc, whichX, T);
        for c = 1:length(conc)
            
            Xtemp = Xfilledall(:,:,c);
            
            SS(c) = sum(sum(Xtemp(indices).^2));
        end
        
        
        f=2*conv;
        iter = 1;
        while iter<max_iter && f>conv
            iter = iter+1;
            SSold = SS; 
            %perform the completion step for each slice  
            for c = 1:length(conc)
                % preprocess = scale and then center
                Xfilled = Xfilledall(:,:,c);
                Xc = Xfilled;
                if center ==1
                    mx = mean(Xc); %columnwise mean 
                    Xc = Xc-ones(m,1)*mx; %centering of the data done - subtract mean from each entry in a column 
                end
                if scale ==1 
                    sj =sqrt(sum((Xc).^2,2));
                    Xc = Xc./(sj*(ones(1,size(Xfilled,2))));
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
                X_pred = St*D';
                if scale ==1
                    X_pred = X_pred.*(sj*(ones(1,size(Xfilled,2))));
                end
                if center ==1
                    X_pred = X_pred+ones(m,1)*mx;
                else 
                    X_pred = St*D';
                end 
                
                
                Xfilled(missing_ind)=X_pred(missing_ind);
                Xfilledall(:,:,c) = Xfilled;
                 
                SS(c) = sum(sum((X_pred(indices)).^2));
            end %end of completion step for concentrations   

            %f = mean(SS);
            
             
            fold = f;
            %check for outliers in each mixture and fix  
            if (mod(iter,55)==0 || fold<f ) && (max_iter-iter)>1000 % every 55 iterations check for outliers 
                for ind = 1:length(missing_ind)
                    tempX = Xfilledall(row(ind),col(ind),:);
                    tempX(abs(tempX)>1e6)=median(tempX);
                    tempScore = (tempX-mean(tempX))./std(tempX);
                    tempX(abs(tempScore)>threshold) = mean(tempX);
                    Xfilledall(row(ind),col(ind),:) = tempX;
                end                  
            end 
            f = mean(abs(SS-SSold)./SSold);
        end
        
        iterations = iter;   
        disp('Iter')
        disp(iter)
        disp('Loss function')
        disp(f) 
        X_predall = Xfilledall;
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