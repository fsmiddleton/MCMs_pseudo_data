function [filenamesave,Xm_boot,Xm_boot2,mse_LOOCV,fns,wmse_LOOCV, X, Xs, conc_interval, filename] = Completion3wayParallelScript()
    % decide on which interval and temperature to evaluate 
    interval = 0.05;
    % comment out the small or all arrays (small contains only alkanes and
    % primary alcohols
    % filename = strcat('HEData3wayPolySmall-',num2str(interval), '-', num2str(T), '.mat');
    % filename = strcat('HEData3wayPolyAll-',num2str(interval), '-', num2str(T), '.mat');
    % load(filename)
    filename = 'HE4wayArrayPolySmall4'; 
    load(filename,  'HE_data_sparse', 'allcomps', 'comps', 'mixtureT','mixture', 'Temps')
    %load('HEData4wayArrayPolySmall-T=4-0.05.mat')
    conc_interval = interval:interval:(1-interval);
    X = HE_data_sparse(:,:,1:length(conc_interval),:);
    Xsign = sign(X);
    Xscale = Xsign.*log(Xsign.*X);
    dim = size(X);
    dim1 = dim(1);
    dim2 = dim(2);
    dim3 = dim(3);
    percmiss = length(find(isnan(X)))/prod(dim)*100;
    percobs = length(find(~isnan(X)))/prod(dim)*100;

    % define all the mixtures in the linear indices of the array, based on the
    % components in the 3-way array file 
    mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 



    % LOOCV to find rank using mse  
    %time the code 

    % declare ranks to be tested 
    fns =[1:2:10];
    maxiter = 10000;
    scale = 1;
    center = 1;
    conv = 1e-10;
    fillmethod = 'avr';
    %orth = 2;
    whichX = 'none';
    %threshold = 0.5*(size(X,3)-1)/sqrt(size(X,3)); 

    %Xscale must have zeros on the diagonal too
    for i =1:dim(4)
        for j = 1:dim(3)
            Xtemp = Xscale(:,:,j,i);
            Xtemp(1:1+size(Xtemp,1):end) = 0; % diagonals are zero
            Xscale(:,:,j,i) = Xtemp;
        end 
    end 

    for Tempind = 1:length(Temps)
        tempX = tril(X(:,:,1,Tempind),-1)+triu(nan(size(X(:,:,1,Tempind))));
        missing_ind{Tempind} = find(~isnan(tempX));
    end 


    %completion algorithm 
    Tempind = 1;%:length(Temps)
    Temperature = Temps(Tempind);
    
    % declare vars for analysis 
    mse_LOOCV = zeros(length(fns),1);
    wmse_LOOCV = zeros(length(fns),1);
    RAD_LOOCV = zeros(length(fns),1); % relative absolute deviation 

     if strcmp(whichX, 'scale') 
        Xs = Xscale;
    else 
        Xs = Xsign;
    end 
    
    
    
    tempX = tril(Xs(:,:,1,Tempind),-1)+triu(nan(size(Xs(:,:,1,Tempind))));
    [row,col] = find(~isnan(tempX));
    filled_ind = find(~isnan(tempX));
    
    error_LOOCV = zeros(length(fns),length(row),dim3);
    error_LOOCV2 = zeros(length(fns),length(row),dim3);
    RAD = zeros(length(fns),length(row),dim3);
    RAD2 = zeros(length(fns),length(row),dim3);
    Xm_boot=zeros( length(fns), length(filled_ind),dim3);
    Xm_boot2=zeros( length(fns), length(filled_ind),dim3);
    %declare vars with size dependent on the array used 
    %loop through ranks
    
    fnind=0;
    for fn = fns 
        disp('fn')
        disp(fn)
        fnloop(1:2) = fn;
        fnloop(3) = fn;
        fnind = fnind + 1; 
        
        parfor k = 1:length(row) %filled_linear_ind must be a row vector for a for loop    
            
            % remove a point from Xs
            X_b = Xs;
            X_b(row(k),col(k),:,Tempind) = nan;
            X_b(col(k),row(k),:,Tempind) = nan;
            if any(find(~isnan(X_b(:,col(k),:,Tempind)))) & any(find(~isnan(X_b(row(k),:,:,Tempind)))) & all(Xs(row(k),col(k),:,Tempind)~=0) %ensure at least one value in each column and row
                disp('k')
                disp(k)
                %[X_pred,F] = missing_parafacThreshold(X_b,fn,maxiter, conv,  orth,conc_interval, Temps, fillmethod,mixtures,whichX,center);
                %Factors{fn,k} = F;
                %[X_pred,iters,F] = missing_parafac4(X_b,fn,100,conv,scale,center,fillmethod,orth, mixtures,conc_interval, whichX,Temps);
                %[X_pred,F] = missing_parafacSimple(X_b,fn,conv,  orth,conc_interval, Temps);  
                [X_pred, iters]=missing_hosvd_par(X_b,fnloop,center,scale,conv,maxiter, fillmethod,mixtures,whichX,conc_interval,Temps);
                
                error_LOOCV(fnind,k, :) = Xs(row(k),col(k),:,Tempind)-X_pred(row(k),col(k),:,Tempind);
                error_LOOCV2(fnind,k, :) = Xs(col(k),row(k),:,Tempind)-X_pred(col(k),row(k),:,Tempind);
                Xm_boot(fnind, k,:) = X_pred(row(k),col(k),:,Tempind);
                Xm_boot2(fnind, k,:) = X_pred(col(k),row(k),:,Tempind);
                if any(Xs(row(k),col(k),:,Tempind)~=0)
                    RAD(fnind,k,:) = error_LOOCV(fnind,k,:)./Xs(row(k),col(k),:,Tempind);
                    RAD2(fnind,k,:) = error_LOOCV2(fnind,k,:)./Xs(col(k),row(k),:,Tempind);
                end
                
            end
        end
        % mse for this composition and rank 
        mse_LOOCV(fnind)= sum(sum(error_LOOCV(fnind,:,:).^2))/length(error_LOOCV(fnind,:,:))+sum(sum(error_LOOCV2(fnind,:,:).^2))/length(error_LOOCV(fnind,:,:));
        wmse_LOOCV(fnind) = find_wmse_error([error_LOOCV(fnind,:,:) error_LOOCV2(fnind,:,:)], length(filled_ind')*2);
        %absolute average deviation
        %RAD_LOOCV(fnind) = sum(sum(abs(RAD(fnind,:,:))))/length(RAD(fnind,:,:))+sum(sum(abs(RAD2(fnind,:,:))))/length(RAD(fnind,:,:));
    end % END FN
    filenamesave = strcat('3wayHOSVD-NoThreshold-centre+scale-Small-concslices-LOOCV-X',whichX,'-maxiter=20000-Temp=',num2str(Temps(Tempind)), '-fillmethod=',fillmethod,'-',  date, '.mat');
    %save(filenamesave) 
 
end 
