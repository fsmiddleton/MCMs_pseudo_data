function [X_pred,iter,F,err] = missing_parafac3(X,fn,max_iter,conv,scale,center, fillmethod)
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
    %
    % Output 
    % S,V,D from X = SVD'
    % St = SV
    % X_pred = filled X with new values for the missing entries
    dim=size(X);
    missing_ind = find(isnan(X));
    filled_ind = find(~isnan(X));
    %choose indices to use for convergence 
    indices=missing_ind;        
   
     
    if any(isnan(X)) % there is missing data
        if size(dim)>3 % more than 3 -way data = 4-way data
            Xfilledini = zeros(dim);
            for i = 1:dim(4)
                Xfilledini(:,:,:,i) = filldata3(X(:,:,:,i),fillmethod);
            end 
        else
            Xfilledini = filldata3(X, fillmethod); 
        end
        f=2*conv;
        iter = 1;
        % PARAFAC options 
         Options(1) = conv;
         Options(2) = 1;
         Options(3) = 0;
         Options(4) = 0;
         Options(5) = 1000;
         Options(6) = max_iter;
         const=[0 0 0]; % orthogonal factors
        %initialise the factors - correct 
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
        % initialises the factors 
        [F,~]=parafac(Xc,fn, Options, const);
        model = nmodel(F);
        % fill predicted values into array
        X_pred = model;
%         X_pred = Xc;
%         X_pred(indices) = model(indices);
        % uncenter and unscale and uncenter 
        if center ==1
            X_pred = X_pred + ones(size(Xfilled,1),1)*mx2;
        end 
        if scale ==1
            X_pred = X_pred.*(sj*(ones(1,size(Xfilled,2))));
        end
        X_pred = X_pred + ones(size(Xfilled,1),1)*mx; 
        
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
%             X_pred = Xc;
%             % fill missing values 
%             X_pred(indices) = model(indices);
            
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
            Xfilledini(indices) = X_pred(indices);
            % calculate sum of squares 
            SS = sum(sum(sum(Xfilledini(indices).^2)));
           
            f = abs(SS-SSold)/(SSold);
            
        end
        err = Xfilled(filled_ind) - X_pred(filled_ind);
    else 
        % no missing data 
        disp('no missing data')

    end %end if  else 
end


function [X_filled, missing] = filldata3(X, method,mixtures,concintervalarray)
    % filldata3 fills a 3-way array with the arithmetic mean of the other
    % values in the array. Used for a 2-way array here
    % Input 
    % X = data array 
    % method = 'avg' - column and row averages 
    %       or 'avc' - column averages 
    %       or 'avr' - row averages 
    %       or 'uni' - unifac predictions 
    %       or 'dia' - upper and lower diagonals used to fill initially 
    % mixtures = dictionary containing components ordered according to their place on the axis 
    % concintervalarray = concentration interval over which to fill data 
    % Output 
    % X_filled = filled data array
    % missing = linear indices of missing data 

    % fill as done for PCA with SVD (averages of each dimension)
    % returns rows and columns with nan elements
    dim = size(X);
    
    if (size(dim,2))==2
        dim(3)=1;
    end 
    missing = isnan(X); %linear indices 

    X_filled = X;
     
    
    if strcmp(method,'avg')
    % f
        % in 2-wy array, do 2-way avg fills
        X_filled = zeros(dim);
        for i = 1:dim(3)
            X_filled(:,:,i) = fill_data(X(:,:,i)); 
        end 
    elseif strcmp(method, 'avc')
        % column averages are used only 
        [m,n]=size(X);
        missing_ind = (isnan(X));
        [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
        X_filled=X;
        X_filled(missing_ind)=0; %fill NaN values with 0
        mean_col = sum(X_filled,1)./(ones(1,n)*m-sum(missing_ind,1));   
        % for all NaN elements that exist, loop through them to replace with means 
        for k =1:length(i) 
            X_filled(i(k),j(k))=mean_col(j(k));
        end  
    elseif strcmp(method, 'avr')
    % row averages are used only 
        [m,n]=size(X);
        missing_ind = (isnan(X));
        [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
        X_filled=X;
        X_filled(missing_ind)=0; %fill NaN values with 0
        mean_row = sum(X_filled,2)./(ones(1,m)*n-sum(missing_ind,2));   
        % for all NaN elements that exist, loop through them to replace with means 
        for k =1:length(i) 
            X_filled(i(k),j(k))=mean_row(i(k));
        end
    elseif strcmp(method, 'dia') % lower and upper diagonals are used in their averages 
        [m,n]=size(X);
        missing_ind = (isnan(X));
        [i, j]=find(isnan(X));% returns rows and columns with nonzero elements
        X_filled=X;
        X_filled(missing_ind)=0; %fill NaN values with 0   
        % for all NaN elements that exist, loop through them to replace with means 
        
        % lower and upper triangular parts of the array
        LX = tril(X_filled,-1);
        UX = triu(X_filled, 1);
        % averages for lower and upper
        mean_rowL = sum(LX,2)./(ones(1,m)*n-sum(missing_ind,2)/2);
        mean_colL = sum(LX,1)./(ones(1,n)*m-sum(missing_ind,1)/2);
        mean_rowU = sum(UX,2)./(ones(1,m)*n-sum(missing_ind,2)/2);
        mean_colU = sum(UX,1)./(ones(1,n)*m-sum(missing_ind,1)/2);
        
        for k =1:length(i) 
            if i(k)>j(k) % lower diagonal
                X_filled(i(k),j(k))=(mean_colL(j(k))+mean_rowL(i(k)))/2;
            else 
                X_filled(i(k),j(k))=(mean_colU(j(k))+mean_rowU(i(k)))/2;
            end 
        end
    else %method ='uni'
        % linearise X to a column vector - entries match mixtures
        X1 = reshape(X(:,:,1),[dim(1)*dim(2),1]); % column vector used for checks and indices
        Xtemp = reshape(X,[dim(1)*dim(2), dim(3)]); % reshaped to column vectors with each column containing a concentration interval
        
        mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
        mix1 = mixtures(:,[1,2]);
        mix2 = mixtures(:,[3,4]);
        % load prediction data and mixtures  
        load('heUNIQUACforT=298.15.mat','he') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
        load('heUNIQUACforT=298.15.mat','mixture') % mixtures for he data - 5151 components
        load('heUNIQUACforT=298.15.mat','conc_interval')
        %convert arrays of concentrations to strings to allow them to be
        %compared 
        conc_unifac = string(conc_interval);
        conc_array = string(concintervalarray);
        conc_array2 = string(1-conc_interval);
        
            
        concindices = find(strcmp(conc_array,conc_unifac));
        concindices2 = find(strcmp(conc_array2,conc_unifac));
        
        %components = possible components in this array, ordered 
        for ind = 1:length(X1)
            if isnan(X1(ind,1))
                % fill with UNIFAC prediction 
                [~,indexpred] = ismember(mixturesarray(ind,:),mixture,'rows');
                
                if indexpred ==0
                    %mixture could be swapped around ie the concentration
                    %is 1-conc as well and is not in the UNIFAC data set 
                    [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixtures,'rows');
                    if indexpred == 0
                        Xtemp(ind,:) = 0;
                        disp(ind)
                    else 
                        Xtemp(ind,:)=he(indexpred,concindices2);
                    end    
                else
                %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                    Xtemp(ind,:)=he(indexpred,concindices);
                end 
            end 
        end 
       
        %reshape X to the 3-way array
        X_filled = reshape(Xtemp,[dim(1), dim(2), dim(3)]);
        %fills remaining few Nan values with averages 
        X_filled = filldata3(X_filled,'avg',mixtures,conc_interval);
    end 
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