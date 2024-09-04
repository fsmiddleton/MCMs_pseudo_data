function [filenamesave,Xm_boot,Xm_boot2,mse_LOOCV,fns,wmse_LOOCV, X, Xs, conc_interval, filename,filled_indices] = Completion3wayParallelScript(fns,fillmethod,threshold,maxiter,filename,fn3)

    load(filename,  'HE_data_sparse', 'allcomps', 'comps', 'mixtureT','mixture', 'Temps')
    interval = 0.05;
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

    mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 

    scale = 1;
    center = 1;
    conv = 1e-10;
    whichX = 'none';

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

    
    % declare vars for analysis 
    mse_LOOCV = zeros(length(fns),1);
    wmse_LOOCV = zeros(length(fns),1);

     if strcmp(whichX, 'scale') 
        Xs = Xscale;
     elseif strcmp(whichX, 'sign')  
        Xs = Xsign;
     else 
         Xs = X;
    end 
    
    tempX = tril(Xs(:,:,1,1),-1)+triu(nan(size(Xs(:,:,1,1))));
    [row,col] = find(~isnan(tempX));
    filled_ind = find(~isnan(tempX));
    
    Xm_boot=zeros(length(Temps),length(fns), length(filled_ind),dim3);
    Xm_boot2=zeros(length(Temps),length(fns), length(filled_ind),dim3);
    
    
    for Tempindex = 1:length(Temps)
        tempX = tril(Xs(:,:,1,Tempindex),-1)+triu(nan(size(Xs(:,:,1,Tempindex))));
        [row,col] = find(~isnan(tempX));
        filled_indices{Tempindex} = find(~isnan(tempX));
        fnind=0;
        for fn = fns 
            fnloop(1:2) = fn;
            fnloop(3) = fn3;
            fnind = fnind + 1; 
            parfor k = 1:length(row) %LOOCV   
                % remove a point from Xs
                X_b = Xs;
                X_b(row(k),col(k),:,Tempindex) = nan;
                X_b(col(k),row(k),:,Tempindex) = nan;
                if any(find(~isnan(X_b(:,col(k),:,Tempindex)))) & any(find(~isnan(X_b(row(k),:,:,Tempindex)))) & all(Xs(row(k),col(k),:,Tempindex)~=0) %ensure at least one value in each column and row
                    [X_pred, iters]=missing_hosvd_par(X_b,fnloop,conv,maxiter, fillmethod,mixtures,whichX,conc_interval,Temps,threshold);
                    Xm_boot(Tempindex,fnind, k,:) = reshape(X_pred(row(k),col(k),:,Tempindex),1,1,[]);
                    Xm_boot2(Tempindex,fnind, k,:) = reshape(X_pred(col(k),row(k),:,Tempindex),1,1,[]);                
                end
            end
        end % END FN
    end 
    filenamesave = strcat('3wayHOSVD-threshold',num2str(threshold),'-fn3=',num2str(fnloop(3)),'-Small-concslices-LOOCV-X',whichX,'-maxiter=20000-fillmethod=',fillmethod,'-',  date, '.mat');
end 
