function [filenamesave,Xm_boot,Xm_boot2,mse_LOOCV,fns,wmse_LOOCV, X, Xs,conc_interval,filename, Factors] = Completion2wayParallelScript()
    % decide on which interval and temperature to evaluate 
    interval = 0.05;
    useprevious = 0;
    
    if useprevious ==1
        load('2waySVD-final-20comps-threshold50par-Xscale-T=298.15-fillmethod=reg-fn=14-18-Nov-2022.mat','X_pred', 'mixtureT', 'filled_ind','row','col')
        
    end
    
    filename = 'HEData3wayPolyMed298.15.mat';
    %filename = strcat('HEData3wayPolySmall-',num2str(interval), '-', num2str(T), '.mat');
    %filename = strcat('HEData3wayPolyAll-',num2str(interval), '-', num2str(T), '.mat');
    load(filename, 'HE_data_sparse',  'comps', 'mixtureT','mixture', 'Temps')
    T = Temps;
    conc_interval = interval:interval:(1-interval);
    
    X = HE_data_sparse;
    Xsign = sign(X);
    Xscale = Xsign.*log(Xsign.*X);
    dim = size(X);
    dim1 = dim(1);
    indend = dim1;
    dim2 = dim(2);
    dim3 = dim(3);
    percmiss = length(find(isnan(X)))/(dim1*dim2*dim3)*100;
    percobs = length(find(~isnan(X)))/(dim1*dim2*dim3)*100;
    
    
    
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
    
    
    tempX = tril(X(:,:,1),-1)+triu(nan(size(X(:,:,1))));
    % fill the X array with previously created answers 
    if useprevious ==1
        %fill Xscale values with those created, redefinition of indend 
        indend = size(X_pred,1);
        Xscale(1:indend,1:indend,:) = X_pred;
        tempX(1:indend,1:indend) = nan;
    end 
    %find the filled values 
    [row,col] = find(~isnan(tempX));
    filled_ind = find(~isnan(tempX));


    % LOOCV to find rank using mse  

    % decl;are ranks to be tested 
    fns =[1:1:10];
    concentrations=conc_interval;

    % vars for missing_parafac3
    maxiter = 40000;
    fillmethod = "avr"; 
    scale = 1;
    center = 1;
    conv = 1e-10;
    whichX = 'none';

    % declare vars for analysis 
    mse_LOOCV = zeros(length(fns),1);
    wmse_LOOCV = zeros(length(fns),1);
    RAD_LOOCV = zeros(length(fns),1); % relative absolute deviation 
    Xm_boot=zeros(length(concentrations), length(fns), length(filled_ind));
    Xm_boot2=zeros(length(concentrations), length(fns), length(filled_ind));

    conc = concentrations;
     if strcmp(whichX, 'scale') 
        Xs = Xscale;
     elseif strcmp(whichX,'sign') 
        Xs = Xsign;
     else
        Xs = X;
    end 

    %declare vars with size dependent on the array used 
    %loop through ranks

    fnind=0;
    for fn = fns 
        disp('fn')
        disp(fn)
        fnind = fnind + 1; 
        error_LOOCV = zeros(length(filled_ind),size(X,3));
        RAD = zeros(length(filled_ind),size(X,3));
        error_LOOCV2 = zeros(length(filled_ind),size(X,3));
        RAD2 = zeros(length(filled_ind),size(X,3));
        parfor k =  1:length(filled_ind) % must be a row vector for a for loop    %56,60,102
            filled_index = filled_ind(k);
            disp('k')
            disp(k)
            % remove a point from Xs
            X_b = Xs;
            X_b(row(k),col(k),:) = nan;
            X_b(col(k),row(k),:) = nan;
            if find(~isnan(X_b(:,col(k),1))) & find(~isnan(X_b(row(k),:,1))) & all(Xs(row(k),col(k),1)~=0) %ensure at least one value in each column and row

                %perform iterative PCA on the slightly more empty matrix 
                [S,V,D,St,X_pred, AllSVs, iterations]=missing_svd_par(X_b,fn,center,scale,conv,maxiter, 1,fillmethod,mixtures,whichX,conc,T);
                %[X_pred,iters,F,err] = missing_svd(X_b,fn,maxiter,conv,scale,center,fillmethod,orth, mixtures,conc, whichX,T);
                %Factors{fn,c,k} = F;
                error_LOOCV(k,:) = Xs(row(k),col(k),:)-X_pred(row(k),col(k),:);
                error_LOOCV2(k,:) = Xs(col(k),row(k),:)-X_pred(col(k),row(k),:);
                Xm_boot(:, fnind, k) = X_pred(row(k),col(k),:);
                Xm_boot2(:, fnind, k) = X_pred(col(k),row(k),:);
                if Xs(row(k),col(k))~=0
                    RAD(k,:) = error_LOOCV(k,:)./reshape(Xs(row(k),col(k),:),1,size(X,3));
                    RAD2(k,:) = error_LOOCV2(k,:)./reshape(Xs(col(k),row(k),:),1,size(X,3));
                end
                Factors{fn,k}{1} = S;
                Factors{fn,k}{2} = V;
                Factors{fn,k}{3} = D;  
            end
        end
        % mse for this composition and rank 
        mse_LOOCV(fnind)= sum(sum(error_LOOCV.^2))/prod(size((find(error_LOOCV))));
        wmse_LOOCV(fnind) = find_wmse_error([error_LOOCV; error_LOOCV2], prod(size((find(error_LOOCV))))*2);
        %absolute average deviation
        RAD_LOOCV(fnind) = sum(sum(abs(RAD)))/prod(size(find(RAD)));

    end % END FN

    filenamesave = strcat('2waySVDNothreshold-',num2str(dim(1)),'comps-addingfrom',num2str(indend),'-threshold60par-LOOCV-X',whichX,'-T=',num2str(T), '-fillmethod=',fillmethod,'-',  date, '.mat');
    %save(filenamesave)
end 
