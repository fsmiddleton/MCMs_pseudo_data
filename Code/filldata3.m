function [X_filled, missing] = filldata3(X, method,mixtures,concintervalarray, whichX, T)
    % filldata3 fills a 3-way array with the arithmetic mean of the other
    % values in the array. Used for a 2-way array here
    % Input 
    % X = data array 
    % method = 'avg' - column and row averages
    %          'avc' - column averages
    %          'avr' - row averages
    %          'row' - linear interpolation in the row
    %          'col' - linear interpolation in the column 
    %          'knn' - linear interpolation in the column and row 
    %          'mix' - nearest 5 of the same type of mixture 
    %          'tri' - triangular filling
    %          'uni' - unifac predictions 
    % mixtures = components ordered according to their place in the linear array (2D) 
    % concinterval = concentration interval over which to fill data 
    % whichX = indicates which form of data it is receiving, necessary for initial guesses: 'scale' or 'sign' or 'none'
    % T = temperature of the system - used for UNIFAC initial guesses 
    %
    % Output 
    % X_filled = filled data array
    % missing = linear indices of missing data 

    dim = size(X);
    
    if (size(dim,2))<3
        dim(3)=1;
        dim(4)=1; 
    elseif size(dim,2)<4
        dim(4) = 1;
    end
    missing = find(isnan(X(:,:,1,1))); %linear indices 
    X_filled = X;
    for Tempind = 1:length(T)
        tempX = tril(X(:,:,1,Tempind),-1)+triu(nan(size(X(:,:,1,Tempind))));
        filled_ind{Tempind} = find(~isnan(tempX));
        [row{Tempind}, col{Tempind}] = find(isnan(tempX));
        missing_ind{Tempind} = find(isnan(tempX));
    end 
    
    if strcmp(method,'avg')
        % in 2-way array, do 2-way avg fills
        X_filled = zeros(dim);
        for j =1:dim(4)
            for i = 1:dim(3)
                X_filled(:,:,i,j) = fill_data(X(:,:,i,j)); 
            end 
        end 
    elseif strcmp(method, 'avc')
        % column averages are used only 
        X_filledtemp = zeros(dim);
        for ind2 =1:dim(4)
            for ind =1:dim(3)
                [m,n]=size(X(:,:,ind,ind2));
                missing_ind = (isnan(X(:,:,ind,ind2)));
                [i, j]=find(isnan(X(:,:,ind)));% returns rows and columns with nonzero elements
                X_filled=X(:,:,ind,ind2);
                X_filled(missing_ind)=0; %fill NaN values with 0
                mean_col = sum(X_filled,1)./(ones(1,n)*m-sum(missing_ind,1));   
                % for all NaN elements that exist, loop through them to replace with means 
                for k =1:length(i) 
                    X_filled(i(k),j(k))=mean_col(j(k));
                end  
               X_filledtemp(:,:,ind,ind2)=X_filled;
            end 
        end 
        X_filled = X_filledtemp;
         
    elseif strcmp(method, 'avr')
    % row averages are used only 
        X_filledtemp = zeros(dim);
        for ind2 =1:dim(4)
            for ind =1:dim(3)
                [m,n]=size(X(:,:,ind,ind2));
                missing_ind = (isnan(X(:,:,ind,ind2)));
                [i, j]=find(isnan(X(:,:,ind,ind2)));% returns rows and columns with nonzero elements
                X_filled=X(:,:,ind,ind2);
                X_filled(missing_ind)=0; %fill NaN values with 0
                mean_row = sum(X_filled,2)./(ones(1,m)*n-sum(missing_ind,2));   
                % for all NaN elements that exist, loop through them to replace with means 
                for k =1:length(i) 
                    X_filled(i(k),j(k))=mean_row(i(k));
                end
                X_filledtemp(:,:,ind,ind2)=X_filled;
            end 
        end 
        X_filled = X_filledtemp;
    elseif strcmp(method,'knn')
        %linear interpolation  
        X_filled = zeros(dim);
        for j =1:dim(4)
            for i = 1:dim(3)
                X_filled1 = fillmissing(X(:,:,i,j), 'Linear', 2, 'EndValues','nearest'); 
                X_filled2 = fillmissing(X(:,:,i,j), 'Linear', 1, 'EndValues','nearest'); 
                %linear interpolation and the end values are set to the
                %nearest non-missing value
                X_filled(:,:,i,j) = (X_filled1 + X_filled2)./2; %average of row and column guesses
            end
        end 
        X_filled = filldata3(X_filled,'avg',mixtures,concintervalarray,whichX,T);
        
    elseif strcmp(method, 'row')
        X_filledtemp = zeros(dim);
        for j = 1:dim(4)
            for ind =1:dim(3)            
                % for all NaN elements that exist, loop through them to replace with means 
                X_filled = X(:,:,ind,j); 
                X_filled = fillmissing(X_filled, 'linear',2,'EndValues','nearest');
                %linear interpolation and the end values are set to the
                %nearest non-missing value
                X_filledtemp(:,:,ind,j)=X_filled;
            end 
        end 
        X_filled = X_filledtemp;
     elseif strcmp(method, 'col')
         %linear interpolation 
        X_filledtemp = zeros(dim);
        for j =1:dim(4)
            for ind =1:dim(3)            
                % for all NaN elements that exist, loop through them to replace with means 
                X_filled = X(:,:,ind,j); 
                X_filled=fillmissing(X_filled, 'linear',1,'EndValues','nearest');
                %linear interpolation and the end values are set to the
                %nearest non-missing value
                X_filledtemp(:,:,ind,j)=X_filled;
            end 
        end 
        X_filled = X_filledtemp;
        
    elseif strcmp(method,'mix') %same mixture type used in the 2-way arrays 
        K=5;
        f1 = mixtures(:,1);
        f2 = mixtures(:,3);
        funcgroups = [f1 f2];
        X_filled = X;
        
        %lower triangular array 
        for ind2 =1:dim(4)
            Xtemp = X;
            [row, col] = find(isnan(tril(Xtemp(:,:,1,ind2),-1)));
            missing = find(isnan(tril(Xtemp(:,:,1,ind2),-1))); %only the missing data on the lower traingular array
            Xtemp = tril(X(:,:,1,ind2),-1);%lower triangular 
            Xtemp(Xtemp==0) = nan; %set zeros to nan to make lower triangular we want
            filled_indlogical = (~isnan(X(:,:,1,ind2))); % finds those entries which need to be filled
            for ind =1:dim(3)            
                % for all NaN elements that exist, loop through them to replace with means
                Xtemp = X(:,:,ind,ind2); % 2-way array 
                %Xtemp = tril(Xtemp,-1);
                for j = 1:length(missing)
                    disp('j')
                    missing_ind = missing(j);
                    mixture = mixtures(missing_ind,:);
                    func1 = mixture(1);
                    func2 = mixture(3);
                    func = [func1 func2];
                    index = (ismember(funcgroups, func,'rows')); %where the same type of mixtures are found in all the mixtures
                    
                    index = reshape(index,dim(1),dim(2));
                    %indices where there are filled values and the same type of
                    %mixture
                    disp(filled_indlogical)
                    index2 = (index.*filled_indlogical);
                    %disp(index2)
                    [row2,col2] = find(index2); %row and column indices of all ones in index2
                    %find those  indices nearest the index of the missing data
                    %point  
                    
                    idx = knnsearch([row2 col2], [row(j) col(j)], 'K', K);
                     if idx
                        for k = 1:length(idx)
                            Xtempk(k) = Xtemp(row2(idx(k)),col2(idx(k)));
                        end 
                        Xtempj = mean(Xtempk);
                        X_filled(row(j),col(j),ind,ind2) = Xtempj;
                     end 
                    
                end 
            end 
        end 
        
        %upper triangular array 
        for ind2 =1:dim(4)
            Xtemp = X;
            [row, col] = find(isnan(triu(Xtemp(:,:,1,ind2),1)));
            missing = find(isnan(triu(Xtemp(:,:,1,ind2),1)));
            Xtemp = triu(X(:,:,1,ind2),1);
            Xtemp(Xtemp==0) = nan;
            filled_indlogical = (~isnan(Xtemp)); % 2-way array 
            
            for ind =1:dim(3)            
                % for all NaN elements that exist, loop through them to replace with means
                Xtemp = X(:,:,ind,ind2);
                Xtemp = triu(Xtemp,1);
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
                    if idx
                        for k = 1:length(idx)
                            Xtempk(k) = Xtemp(row2(idx(k)),col2(idx(k)));
                        end 
                        Xtempj = mean(Xtempk);
                        X_filled(row(j),col(j),ind,ind2) = Xtempj;
                    end 
                end 
            end 
        end
        X_filled = filldata3(X_filled,'avg',mixtures,concintervalarray,whichX,T);
        
    
    elseif strcmp(method, 'tri') % lower and upper diagonals are used in their averages
        X_filledtemp = zeros(dim);
        for ind2=1:dim(4)
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
        end 
        X_filled = X_filledtemp;
        
    else %method ='uni'
        
        if dim(4)~=1 % 4 dimensions - 4th dimension changed to 1
            for count = 1:length(T)
                missing_indices = missing_ind{count};
                rows = row{count};
                cols = col{count};
                % linearise X to a column vector - entries match mixtures
                X1 = reshape(X(:,:,1,count),[dim(1),[]]); % column vector used for checks and indices
                Xtemp = (X(:,:,:,count)); % reshaped to column vectors with each column containing a concentration interval

                mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
                mix1 = mixtures(:,[1,2]);
                mix2 = mixtures(:,[3,4]);
                % load prediction data and mixtures  
                load(strcat('heUNIFACforT=',num2str(T(count)),'.mat'),'he','mixture','conc_interval') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
                %convert arrays of concentrations to strings to allow them to be
                %compared 
                conc_unifac = string(conc_interval);
                conc_array =string(concintervalarray);
                mixture = mixture'; % rows can now be compared 

                concindices2 = (intersect(conc_unifac,conc_array));
                
                [~,concindices2] = ismember(concindices2,conc_unifac);
                
                %components = possible components in this array, ordered 
                
                for ind = 1:length(missing_indices)
                    if isnan(X1(ind,1)) 
                        % fill with UNIFAC prediction 
                        %disp(ind)
                        index = missing_indices(ind);
                        [~,indexpred] = ismember(mixturesarray(index,:),mixture,'rows');

                        if indexpred ==0
                            %mixture could be swapped around ie the concentration
                            %is 1-conc as well and is not in the UNIFAC data set 
                            [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                            if indexpred == 0
                                Xtemp(rows(ind), cols(ind),:) = 0;
                                Xtemp(cols(ind), rows(ind),:) = 0;
                            else 
                                if indexpred> size(he,2)
                                    [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                                end 
                                temp = he(concindices2, indexpred);

                                if strcmp(whichX, 'scale')
                                    Xtemp(rows(ind), cols(ind),:)=sign(temp).*log(sign(temp).*temp);
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)=sign(temp).*log(sign(temp).*temp);
                                elseif strcmp(whichX, 'sign') 
                                    Xtemp(rows(ind), cols(ind),:)= sign(temp);
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)= sign(temp);
                                else 
                                    Xtemp(rows(ind), cols(ind),:)=temp;
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)=temp;
                                end 
                            end    
                        else
                        %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                            if indexpred> size(he,2)
                                [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                            end     

                            temp = he(concindices2, indexpred);

                            if strcmp(whichX, 'scale')
                                Xtemp(rows(ind), cols(ind),:)=sign(temp).*log(sign(temp).*temp);
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)=sign(temp).*log(sign(temp).*temp);
                            elseif strcmp(whichX, 'sign') 
                                Xtemp(rows(ind), cols(ind),:)= sign(temp);
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)= sign(temp);
                            else 
                                Xtemp(rows(ind), cols(ind),:)=temp;
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)=temp;
                            end  
                        end 
                    end 
                end 
                %reshape X to the 3-way array
                X_filledtemp = reshape(Xtemp,[dim(1), dim(2), dim(3)]);
                %fills remaining few Nan values with averages 
                X_filled(:,:,:,count) = filldata3(X_filledtemp,'avg',mixtures,conc_interval,whichX,T);
            end 
        elseif dim(3)==1 %2 dimensions 
                count =1;
                missing_indices = missing_ind{count};
                rows = row{count};
                cols = col{count};
                X1 = reshape(X(:,:,1,count),[dim(1),[]]); % column vector used for checks and indices
                Xtemp = (X(:,:,:,count)); % reshaped to column vectors with each column containing a concentration interval

                mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
                mix1 = mixtures(:,[1,2]);
                mix2 = mixtures(:,[3,4]);
                % load prediction data and mixtures  
                load(strcat('heUNIFACforT=',num2str(T(count)),'.mat'),'he','mixture','conc_interval') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
                %convert arrays of concentrations to strings to allow them to be
                %compared 
                conc_unifac = string(conc_interval);
                conc_array =string(concintervalarray);
                mixture = mixture'; % rows can now be compared 

                concindices2 = (intersect(conc_unifac,conc_array));
                
                [~,concindices2] = ismember(concindices2,conc_unifac);
                
                %components = possible components in this array, ordered 
                
                for ind = 1:length(missing_indices)
                    if isnan(X1(ind,1)) 
                        % fill with UNIFAC prediction 
                        %disp(ind)
                        index = missing_indices(ind);
                        [~,indexpred] = ismember(mixturesarray(index,:),mixture,'rows');

                        if indexpred ==0
                            %mixture could be swapped around ie the concentration
                            %is 1-conc as well and is not in the UNIFAC data set 
                            [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                            if indexpred == 0
                                Xtemp(rows(ind), cols(ind),:) = 0;
                                Xtemp(cols(ind), rows(ind),:) = 0;
                            else 
                                if indexpred> size(he,2)
                                    [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                                end 
                                temp = he(concindices2, indexpred);

                                if strcmp(whichX, 'scale')
                                    Xtemp(rows(ind), cols(ind),:)=sign(temp).*log(sign(temp).*temp);
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)=sign(temp).*log(sign(temp).*temp);
                                elseif strcmp(whichX, 'sign') 
                                    Xtemp(rows(ind), cols(ind),:)= sign(temp);
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)= sign(temp);
                                else 
                                    Xtemp(rows(ind), cols(ind),:)=temp;
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)=temp;
                                end 
                            end    
                        else
                        %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                            if indexpred> size(he,2)
                                [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                            end     

                            temp = he(concindices2, indexpred);

                            if strcmp(whichX, 'scale')
                                Xtemp(rows(ind), cols(ind),:)=sign(temp).*log(sign(temp).*temp);
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)=sign(temp).*log(sign(temp).*temp);
                            elseif strcmp(whichX, 'sign') 
                                Xtemp(rows(ind), cols(ind),:)= sign(temp);
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)= sign(temp);
                            else 
                                Xtemp(rows(ind), cols(ind),:)=temp;
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)=temp;
                            end  
                        end 
                    end 
                end 
                %reshape X to the 3-way array
                X_filledtemp = reshape(Xtemp,[dim(1), dim(2)]);
                %fills remaining few Nan values with averages 
                X_filled = filldata3(X_filledtemp,'avg',mixtures,conc_interval,whichX,T);
        else %3 dimensions
            % linearise X to a column vector - entries match mixtures
                count =1;
                missing_indices = missing_ind{count};
                rows = row{count};
                cols = col{count};
                X1 = reshape(X(:,:,1,count),dim(1),[]); % column vector used for checks and indices
                Xtemp = (X(:,:,:,count)); % reshaped to column vectors with each column containing a concentration interval

                mixturesarray = mixtures; % each row contains the mixture of that X (if X were a 2-way array)
                mix1 = mixtures(:,[1,2]);
                mix2 = mixtures(:,[3,4]);
                % load prediction data and mixtures  
                load(strcat('heUNIFACforT=',num2str(T(count)),'.mat'),'he','mixture','conc_interval') % he in J/mol for composition of 0.01 to 0.99 for 5151 components
                %convert arrays of concentrations to strings to allow them to be
                %compared 
                conc_unifac = string(conc_interval);
                conc_array =string(concintervalarray);
                mixture = mixture'; % rows can now be compared 

                concindices2 = (intersect(conc_unifac,conc_array));
                
                [~,concindices2] = ismember(concindices2,conc_unifac);
                
                %components = possible components in this array, ordered 
                
                for ind = 1:length(missing_indices)
                    if isnan(X1(rows(ind), cols(ind),1)) 
                        % fill with UNIFAC prediction 
                        %disp(ind)
                        index = missing_indices(ind);
                        [~,indexpred] = ismember(mixturesarray(index,:),mixture,'rows');

                        if indexpred ==0
                            %mixture could be swapped around ie the concentration
                            %is 1-conc as well and is not in the UNIFAC data set 
                            [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                            if indexpred == 0
                                Xtemp(rows(ind), cols(ind),:) = 0;
                                Xtemp(cols(ind), rows(ind),:) = 0;
                            else 
                                if indexpred> size(he,2)
                                    [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                                end 
                                temp = he(concindices2, indexpred);

                                if strcmp(whichX, 'scale')
                                    Xtemp(rows(ind), cols(ind),:)=sign(temp).*log(sign(temp).*temp);
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)=sign(temp).*log(sign(temp).*temp);
                                elseif strcmp(whichX, 'sign') 
                                    Xtemp(rows(ind), cols(ind),:)= sign(temp);
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)= sign(temp);
                                else 
                                    Xtemp(rows(ind), cols(ind),:)=temp;
                                    temp = flip(temp);
                                    Xtemp(cols(ind), rows(ind),:)=temp;
                                end 
                            end    
                        else
                        %indexpred = find(ismember(mixturestemp(ind),mixture, 'rows'));%find mixture that is missing in the unifac prediction mixtures                
                            if indexpred> size(he,2)
                                [~,indexpred] = ismember([mix2(index,:) mix1(index,:)],mixture,'rows');
                            end     

                            temp = he(concindices2, indexpred);

                            if strcmp(whichX, 'scale')
                                Xtemp(rows(ind), cols(ind),:)=sign(temp).*log(sign(temp).*temp);
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)=sign(temp).*log(sign(temp).*temp);
                            elseif strcmp(whichX, 'sign') 
                                Xtemp(rows(ind), cols(ind),:)= sign(temp);
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)= sign(temp);
                            else 
                                Xtemp(rows(ind), cols(ind),:)=temp;
                                temp = flip(temp);
                                Xtemp(cols(ind), rows(ind),:)=temp;
                            end  
                        end 
                    end 
                end 
                %reshape X to the 3-way array
                X_filledtemp = reshape(Xtemp,[dim(1), dim(2), dim(3)]);
                %fills remaining few Nan values with averages 
                X_filled = filldata3(X_filledtemp,'avg',mixtures,conc_interval,whichX,T);
        end 
        
    end
    
    %check for any missing data and fill with averages if they exist 
    if any(isnan(X))
        if length(dim)>3
            for j=1:dim(4)
                for i=1:dim(3)
                    X_filled(:,:,i,j)=fill_data(X_filled(:,:,i,j));
                end 
            end 
        else 
            for i=1:dim(3)
                X_filled(:,:,i)=fill_data(X_filled(:,:,i));
            end 
        end 
        
    end
    X_filled = reshape(X_filled, size(X));
    
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
    %remove nan vals 
    mean_row(isnan(mean_row)) =0;
    mean_col(isnan(mean_col)) =0;
    % for all NaN elements that exist, loop through them to replace with means 
    for k =1:length(i) 
        X_filled(i(k),j(k))=(mean_row(i(k))+mean_col(j(k)))/2;
    end  
end 
