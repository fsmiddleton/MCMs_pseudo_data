function [X_filled, missing] = filldata3(X, method,mixtures,concintervalarray, whichX, T)
    % filldata3 fills a 3-way array with the arithmetic mean of the other
    % values in the array. Used for a 2-way array here
    % Input 
    % X = data array 
    % method = 'avg' - column and row averages
    %          'dia' - triangular filling
    %          'avc' - column averages
    %          'avr' - row averages
    %          'reg' - nearest values in the row
    %          'rec' - nearest values in the column 
    %          'reb' - nearest values in the column and row 
    %          'mix' - nearest values of the same type of mixture 
    %       or 'uni' - unifac predictions 
    % mixtures = components ordered according to their place in the linear array (2D) 
    % concinterval = concentration interval over which to fill data 
    % whichX = indicates which form of data it is receiving, necessary for initial guesses: 'scale' or 'sign' or 'none'
    % T = temperature of the system - used for UNIFAC initial guesses 
    %
    % Output 
    % X_filled = filled data array
    % missing = linear indices of missing data 

    % fill as done for PCA with SVD (averages of each dimension)
    % returns rows and columns with nan elements
    dim = size(X);
    
    if (size(dim,2))==2
        dim(3)=1;
    end 
    missing = find(isnan(X(:,:,1))); %linear indices 
    

    X_filled = X;
     
    
    if strcmp(method,'avg')
    % f
        % in 2-way array, do 2-way avg fills
        X_filled = zeros(dim);
        for i = 1:dim(3)
            X_filled(:,:,i) = fill_data(X(:,:,i)); 
        end 
    elseif strcmp(method,'reb')
        X_filled = zeros(dim);
        for i = 1:dim(3)
            X_filled1 = fillmissing(X(:,:,i), 'Linear', 2, 'EndValues','nearest'); 
            X_filled2 = fillmissing(X(:,:,i), 'Linear', 1, 'EndValues','nearest'); 
            X_filled(:,:,i) = (X_filled1 + X_filled2)./2;
        end
    elseif strcmp(method, 'reg')
        X_filledtemp = zeros(dim);
        for ind =1:dim(3)            
            % for all NaN elements that exist, loop through them to replace with means 
            X_filled = X(:,:,ind); 
            X_filled=fillmissing(X_filled, 'linear',2,'EndValues','nearest');
             
           X_filledtemp(:,:,ind)=X_filled;
        end 
        X_filled = X_filledtemp;
     elseif strcmp(method, 'rec')
        X_filledtemp = zeros(dim);
        for ind =1:dim(3)            
            % for all NaN elements that exist, loop through them to replace with means 
            X_filled = X(:,:,ind); 
            X_filled=fillmissing(X_filled, 'linear',1,'EndValues','nearest');
             
           X_filledtemp(:,:,ind)=X_filled;
        end 
        X_filled = X_filledtemp;
        
    elseif strcmp(method,'mix') %same mixture type used in the 2-way arrays 
        K=3;
        [row, col] = find(isnan(X(:,:,1)));
        Xtemp = X;
        Xtemp(Xtemp==0) = nan;
        filled_indlogical = (~isnan(Xtemp(:,:,1))); % 2-way array - 1 = not NaN
        
        f1 = mixtures(:,1);
        f2 = mixtures(:,3);
        funcgroups = [f1 f2];
        X_filled = X;
        for ind =1%:dim(3)            
            % for all NaN elements that exist, loop through them to replace with means
            Xtemp = X(:,:,ind);
            for j = 1:length(missing)
                %disp('j')
                missing_ind = missing(j);
                mixture = mixtures(missing_ind,:);
                func1 = mixture(1);
                func2 = mixture(3);
                func = [func1 func2];
                index = (ismember(funcgroups, func,'rows')); %where the same type of mixtures are found in all the mixtures
                index = reshape(index,dim(1),dim(2));
                %indices where there are filled values and the same type of
                %mixture
                index2 = (index.*filled_indlogical);
                [row2,col2] = find(index2); %row and column indices of all ones in index2
                %find those  indices nearest the index of the missing data
                %point                 
                idx = knnsearch([row2 col2], [row(j) col(j)], 'K', K);

                for k = 1:K
                    Xtempk(k) = Xtemp(row2(idx(k)),col2(idx(k)));
                end 
                Xtempj = mean(Xtempk);
                X_filled(row(j),col(j),ind) = Xtempj;
            end 
        end 
    elseif strcmp(method, 'avc')
        % column averages are used only 
        X_filledtemp = zeros(dim);
        for ind =1:dim(3)
            [m,n]=size(X(:,:,ind));
            missing_ind = (isnan(X(:,:,ind)));
            [i, j]=find(isnan(X(:,:,ind)));% returns rows and columns with nonzero elements
            X_filled=X(:,:,ind);
            X_filled(missing_ind)=0; %fill NaN values with 0
            mean_col = sum(X_filled,1)./(ones(1,n)*m-sum(missing_ind,1));   
            % for all NaN elements that exist, loop through them to replace with means 
            for k =1:length(i) 
                X_filled(i(k),j(k))=mean_col(j(k));
            end  
           X_filledtemp(:,:,ind)=X_filled;
        end 
        X_filled = X_filledtemp;
         
    elseif strcmp(method, 'avr')
    % row averages are used only 
        X_filledtemp = zeros(dim);
        for ind =1:dim(3)
            [m,n]=size(X(:,:,ind));
            missing_ind = (isnan(X(:,:,ind)));
            [i, j]=find(isnan(X(:,:,ind)));% returns rows and columns with nonzero elements
            X_filled=X(:,:,ind);
            X_filled(missing_ind)=0; %fill NaN values with 0
            mean_row = sum(X_filled,2)./(ones(1,m)*n-sum(missing_ind,2));   
            % for all NaN elements that exist, loop through them to replace with means 
            for k =1:length(i) 
                X_filled(i(k),j(k))=mean_row(i(k));
            end
            X_filledtemp(:,:,ind)=X_filled;
        end 
        X_filled = X_filledtemp;
         
    elseif strcmp(method, 'dia') % lower and upper diagonals are used in their averages
        X_filledtemp = zeros(dim);
        for ind = 1:dim(3)
            [m,n]=size(X(:,:,ind));
            missing_ind = (isnan(X(:,:,ind)));
            [i, j]=find(isnan(X(:,:,ind)));% returns rows and columns with nonzero elements
            X_filled=X(:,:,ind);
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
            X_filledtemp(:,:,ind)=X_filled;
        end 
        X_filled = X_filledtemp;
        
    else %method ='uni'
        
        if size(dim,2)>3
            for count = 1:length(T)
                % linearise X to a column vector - entries match mixtures
                X1 = reshape(X(:,:,1),[dim(1)*dim(2),1]); % column vector used for checks and indices
                Xtemp = reshape(X(:,:,:,count),[dim(1)*dim(2), dim(3)]); % reshaped to column vectors with each column containing a concentration interval

                mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
                mix1 = mixtures(:,[1,2]);
                mix2 = mixtures(:,[3,4]);
                % load prediction data and mixtures  
                load(strcat('heUNIFACforT=',num2str(T(count)),'.mat'),'he') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
                load(strcat('heUNIFACforT=',num2str(T(count)),'.mat'),'mixture') % mixtures for he data - 5151 components
                load(strcat('heUNIFACforT=',num2str(T(count)),'.mat'),'conc_interval')
                %convert arrays of concentrations to strings to allow them to be
                %compared 
                conc_unifac = string(conc_interval);
                conc_array =string(concintervalarray);
                conc_unifac2 = (1-conc_interval);
                mixture = mixture'; % rows can now be compared 

                concindices2 = (intersect(conc_unifac,conc_array));
                
                [~,concindices2] = ismember(concindices2,conc_unifac);
                
                %components = possible components in this array, ordered 
                
                for ind = 1:length(X1)
                    if isnan(X1(ind,1)) 
                        % fill with UNIFAC prediction 
                        %disp(ind)
                        [~,indexpred] = ismember(mixturesarray(ind,:),mixture,'rows');

                        if indexpred ==0
                            %mixture could be swapped around ie the concentration
                            %is 1-conc as well and is not in the UNIFAC data set 
                            [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                            if indexpred == 0
                                Xtemp(ind,:) = 0;
                                
                            else 
                                if indexpred> size(he,2)
                                    [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                                end 
                                temp = he(concindices2, indexpred);

                                if strcmp(whichX, 'scale')
                                    Xtemp(ind,:)=sign(temp).*log(sign(temp).*temp);
                                elseif strcmp(whichX, 'sign') 
                                    Xtemp(ind,:)= sign(temp);
                                else 
                                    Xtemp(ind,:)=temp;
                                end 
                            end    
                        else
                        %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                            if indexpred> size(he,2)
                                [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                            end     

                            temp = he(concindices2, indexpred);

                            if strcmp(whichX, 'scale')
                                Xtemp(ind,:)=sign(temp).*log(sign(temp).*temp);
                            elseif strcmp(whichX, 'sign') 
                                Xtemp(ind,:)= sign(temp);
                            else
                                Xtemp(ind,:)=temp;
                            end 
                        end 
                    end 
                end 
                %reshape X to the 3-way array
                X_filledtemp = reshape(Xtemp,[dim(1), dim(2), dim(3)]);
                %fills remaining few Nan values with averages 
                X_filled(:,:,:,count) = filldata3(X_filledtemp,'avg',mixtures,conc_interval);
            end 
        elseif size(dim)==2
                Xtemp = reshape(X(:,:),[dim(1)*dim(2),1]);
                X1 = reshape(X(:,:,1),[dim(1)*dim(2),1]); % column vector used for checks and indices
                Xtemp = reshape(X,[dim(1)*dim(2), dim(3)]); % reshaped to column vectors with each column containing a concentration interval

                mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
                mix1 = mixtures(:,[1,2]);
                mix2 = mixtures(:,[3,4]);
                % load prediction data and mixtures  
                load(strcat('heUNIFACforT=',num2str(T),'.mat'),'he') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
                load(strcat('heUNIFACforT=',num2str(T),'.mat'),'mixture') % mixtures for he data - 5151 components
                load(strcat('heUNIFACforT=',num2str(T),'.mat'),'conc_interval')
                %convert arrays of concentrations to strings to allow them to be
                %compared 
                conc_unifac = string(conc_interval);
                conc_array =string(concintervalarray);
                
                mixture = mixture'; % rows can now be compared 

                concindices2 = (intersect(conc_unifac,conc_array));
                
                [~,concindices2] = ismember(concindices2,conc_unifac);
                
                %components = possible components in this array, ordered 
                
                for ind = 1:length(X1)
                    if isnan(X1(ind,1)) 
                        % fill with UNIFAC prediction 
                        %disp(ind)
                        [~,indexpred] = ismember(mixturesarray(ind,:),mixture,'rows');

                        if indexpred ==0
                            %mixture could be swapped around ie the concentration
                            %is 1-conc as well and is not in the UNIFAC data set 
                            [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                            if indexpred == 0
                                Xtemp(ind,:) = 0;
                                
                            else 
                                if indexpred> size(he,2)
                                    [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                                end 
                                temp = he(concindices2, indexpred);

                                if strcmp(whichX, 'scale')
                                    Xtemp(ind,:)=sign(temp).*log(sign(temp).*temp);
                                elseif strcmp(whichX, 'sign') 
                                    Xtemp(ind,:)= sign(temp);
                                else 
                                    Xtemp(ind,:)=temp;
                                end 
                            end    
                        else
                        %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                            if indexpred> size(he,2)
                                [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                            end     

                            temp = he(concindices2, indexpred);

                            if strcmp(whichX, 'scale')
                                Xtemp(ind,:)=(sign(temp).*log(sign(temp).*temp))';
                            elseif strcmp(whichX, 'sign') 
                                Xtemp(ind,:)= sign(temp);
                            else
                                Xtemp(ind,:)=temp;
                            end 
                        end 
                    end 
                end 
                %reshape X to the 3-way array
                X_filledtemp = reshape(Xtemp,[dim(1), dim(2)]);
                %fills remaining few Nan values with averages 
                X_filled = filldata3(X_filledtemp,'avg',mixtures,conc_interval);
        else 
            % linearise X to a column vector - entries match mixtures
                X1 = reshape(X(:,:,1),[dim(1)*dim(2),1]); % column vector used for checks and indices
                Xtemp = reshape(X,[dim(1)*dim(2), dim(3)]); % reshaped to column vectors with each column containing a concentration interval

                mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
                mix1 = mixtures(:,[1,2]);
                mix2 = mixtures(:,[3,4]);
                % load prediction data and mixtures  
                load(strcat('heUNIFACforT=',num2str(T),'.mat'),'he') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
                load(strcat('heUNIFACforT=',num2str(T),'.mat'),'mixture') % mixtures for he data - 5151 components
                load(strcat('heUNIFACforT=',num2str(T),'.mat'),'conc_interval')
                %convert arrays of concentrations to strings to allow them to be
                %compared 
                conc_unifac = string(conc_interval);
                %disp(concintervalarray)
                conc_array =string(concintervalarray);
                conc_unifac2 = (1-conc_interval);
                mixture = mixture'; % rows can now be compared 
                
                concindices2 = (intersect(conc_unifac,conc_array));
                %disp(concindices2);
                [~,concindices2] = ismember(concindices2,conc_unifac);
                %disp(concindices2);
                %components = possible components in this array, ordered 
                
                for ind = 1:length(X1)
                    if isnan(X1(ind,1)) 
                        % fill with UNIFAC prediction 
                        [~,indexpred] = ismember(mixturesarray(ind,:),mixture,'rows');

                        if indexpred ==0
                            %mixture could be swapped around ie the concentration
                            %is 1-conc as well and is not in the UNIFAC data set 
                            [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                            if indexpred == 0
                                Xtemp(ind,:) = 0;
                                
                            else 
                                if indexpred> size(he,2)
                                    [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                                end 
                                temp = he(concindices2, indexpred);

                                if strcmp(whichX, 'scale')
                                    Xtemp(ind,:)=sign(temp).*log(sign(temp).*temp);
                                elseif strcmp(whichX, 'sign') 
                                    Xtemp(ind,:)= sign(temp);
                                else 
                                    Xtemp(ind,:)=temp;
                                end 
                            end    
                        else
                        %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                            if indexpred> size(he,2)
                                [~,indexpred] = ismember([mix2(ind,:) mix1(ind,:)],mixture,'rows');
                            end     

                            temp = he(concindices2, indexpred);

                            if strcmp(whichX, 'scale')
                                Xtemp(ind,:)=(sign(temp).*log(sign(temp).*temp))';
                            elseif strcmp(whichX, 'sign') 
                                Xtemp(ind,:)= sign(temp);
                            else
                                Xtemp(ind,:)=temp;
                            end 
                        end 
                    end 
                end 
                %reshape X to the 3-way array
                X_filledtemp = reshape(Xtemp,[dim(1), dim(2), dim(3)]);
                %fills remaining few Nan values with averages 
                X_filled = filldata3(X_filledtemp,'avg',mixtures,conc_interval);
        end 
        
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
