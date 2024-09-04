function [X_predall, iterations]=missing_hosvd_par(X,fn,conv,maxiter, fillmethod,mixtures,whichX,conc,Temps,thresh)
    % Fill a matrix of missing data using PCA with SVD and a given number of
    % PCs. Can also handle non-missing data. Missing data is handled as NaN
    % values 
    %
    % Input 
    % X = 4-way array 
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
    
    dim = size(X);
    indices = find(isnan(reshape(X(:,:,1,:),dim(1),dim(2),dim(4))));
    broprocess = 1;
    threshold = thresh/100*(size(X,3)*2-1)/sqrt(size(X,3)*2);
    cent = [1 1 1];
    scale = [1 1 1];
    %choose indices to use for convergence 
    SS = zeros(length(conc),1); 
    X(X==0)=nan;
    if any(isnan(X(:,:,1,:))) % there is missing data 
        Xfilledall = zeros(dim);
        for i = 1:dim(4)
            Xfilledall(:,:,:,i) = filldata3(X(:,:,:,i),fillmethod, mixtures,conc, whichX, Temps(i));
        end
        
        %find first SS using initial guesses
        for c = 1:length(conc)
            for Tind = length(Temps)
                Xtemp = Xfilledall(:,:,c,Tind);
                tempX = X(:,:,c,Tind);
                missing_ind{Tind} = find(isnan(tempX));
                SS(Tind,c) = (sum(sum(Xtemp(missing_ind{Tind}).^2)));
            end 
        end
        SS = mean(mean(SS));
        
        f=2*conv;
        iter = 1;
        while iter<maxiter && f>conv
            iter = iter+1;
            SSold = SS; 
            %perform the completion step for each slice  
            for c = 1:length(conc)
                %extract correct 3-way array 
                Xloop = Xfilledall(:,:,c,:);
                A = reshape(Xloop,dim(1),dim(2),dim(4));%3way array 
                dimA = size(A);
                %clause to not allow inf or nan values to continue 
                if any(isnan(A),'all') || any(isinf(A),'all')
                    A(isinf(A))=nan;
                    Atemp = filldata3(reshape(A,dim(1),dim(2),1,dim(4)),'avg',mixtures,c, whichX, Temps);
                    A=reshape(Atemp,dim(1),dim(2),dim(4));
                end
                 %A1 prediction 
                    if broprocess 
                        [Ac,mX,sX]=nprocess(A,cent,scale);
                    else 
                        [A1, ~,~]=nshape(A);
                        if center(1) ==1 
                            mx = mean(A1);
                            A1 = A1-ones(size(A1,1),1)*mx; %columnwise mean 
                        end 
                        if scale(1) ==1 
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
                        [A1_pred]=nprocess(A1_pred,cent,scale,mX,sX,-1);
                    else 
                        if scale(1) ==1
                            A1_pred = A1_pred.*(sj*ones(1,size(A1,2)));
                        end
                        if center(1) ==1
                            A1_pred = A1_pred+ones(size(A1,1),1)*mx;
                        end 
                        %refold
                        A1_pred = reshape(A1_pred,dimA);
                    end 
                    X_pred = A1_pred;
                

                % fill missing values with predictions 
                Xloop(indices) = X_pred(indices);
                Xfilledall(:,:,c,:) = Xloop;
               
            end %end of conc_interval  
            
            
            fold = f;
            %finding the convergence criterion 
            for c = 1:length(conc)
                for Tind = length(Temps)
                    Xtemp = Xfilledall(:,:,c,Tind);
                    SS(Tind,c) = (sum(sum(Xtemp(missing_ind{Tind}).^2)));
                end 
            end
            SS = mean(mean(SS));
            f = abs((SS-SSold)/SSold);
            
            %thresholding - evaluation 
            if (mod(iter,55)==0 || fold<f )&& (maxiter-iter)>1000 %every 55 iterations, or if the convergence criterion became worse 
                for Tempindloop = 1:length(Temps)
                    %tempX = X(:,:,1,Tempindloop);
                    tempX = tril(X(:,:,1,Tempindloop),-1)+triu(zeros(size(X(:,:,1,Tempindloop))));
                    [rowloop,colloop] = find(isnan(tempX)); %find missing values in this Temperature slice   
                    for ind = 1:length(rowloop)
                        %extract predictions for this missing value 
                        tempX = [Xfilledall(rowloop(ind),colloop(ind),:,Tempindloop); Xfilledall(colloop(ind),rowloop(ind),:,Tempindloop)];
                        tempX = reshape(tempX,numel(tempX),[]); % reshape to a column vector
                        if any(abs(tempX)>5e5)  %don't fix on the last iteration 
                            tempX(abs(tempX)>5e5)=median(tempX); %set great outliers to nan to recalculate 
                            tempScore = (tempX-mean(tempX))./std(tempX); %calculate scores
                            tempX(abs(tempScore)>threshold) = median(tempX); %remove outliers 
                            Xfilledall(rowloop(ind),colloop(ind),:,Tempindloop) = tempX(1:length(conc));
                            Xfilledall(colloop(ind),rowloop(ind),:,Tempindloop) = tempX(length(conc)+1:length(conc)*2);
                        end 
                    end 
                end 
            end 
        
            
        end
            iterations = iter;   
            X_predall = Xfilledall;
    else 
        disp('No Missing Data')
    end %end if  else 
end % end function