function [filenamesave,Xm_boot,Xm_boot2,fns, X, Xs,conc_interval,filename, filled_ind] = completion_2way_par(fns,fillmethod,maxiter,filename,thresholdperc)

    interval = 0.05;

    load(filename, 'HE_data_sparse',  'comps', 'mixtureT','mixture', 'Tlower')
    T = Tlower+1;
    conc_interval = interval:interval:(1-interval);
    
    X = HE_data_sparse;
    dim = size(X);
    dim1 = dim(1);
    indend = dim1;
    dim2 = dim(2);
    dim3 = dim(3);
    percmiss = length(find(isnan(X)))/(dim1*dim2*dim3)*100;
    percobs = length(find(~isnan(X)))/(dim1*dim2*dim3)*100;
    
    mixtures = zeros(size(comps,1)^2,4);
    index = 0;
    for i = 1:length(comps)
        for j = 1:length(comps)
            index = index+1;
            mixtures(index,:) = [comps(i,:) comps(j,:)];
        end
    end 
    
    tempX = tril(X(:,:,1),-1)+triu(nan(size(X(:,:,1))));
    [row,col] = find(~isnan(tempX));
    filled_ind = find(~isnan(tempX));
    concentrations=conc_interval;

    scale = 1;
    center = 1;
    conv = 1e-10;

    Xm_boot=zeros(length(concentrations), length(fns), length(filled_ind));
    Xm_boot2=zeros(length(concentrations), length(fns), length(filled_ind));

    conc = concentrations;
    
    Xs = X;
    
    %loop through ranks
    fnind=0;
    for fn = fns 
        disp('fn')
        disp(fn)
        fnind = fnind + 1; 
        
        parfor k =  1:length(filled_ind) % LOOCV
            % remove a point from Xs
            X_b = Xs;
            X_b(row(k),col(k),:) = nan;
            X_b(col(k),row(k),:) = nan;
            if find(~isnan(X_b(:,col(k),1))) & find(~isnan(X_b(row(k),:,1))) & all(Xs(row(k),col(k),1)~=0) %ensure at least one value in each column and row
                %perform iterative PCA on the slightly more empty matrix 
                [S,V,D,St,X_pred, AllSVs]=missing_svd_par(X_b,fn,center,scale,conv,maxiter, 1,fillmethod,mixtures,conc,T, thresholdperc);
                Xm_boot(:, fnind, k) = X_pred(row(k),col(k),:);
                Xm_boot2(:, fnind, k) = X_pred(col(k),row(k),:);
            end
        end
    end % end fn
    filenamesave = strcat('2waySVD-',num2str(dim(1)),'comps-threshold',num2str(thresholdperc), 'par-LOOCV-T=',num2str(T), '-fillmethod=',fillmethod,'-',  date, '.mat');
end 
