function  [fn,msefill,wmsefill,RAD, X_pred,cutoff,SVs,iterations] = solveGavish(Xs, dim1, dim2, conv, iter, mixtures, conc,T)
    % Gavish Hard thresholding used to find the rank of the matrix using PCA
    % with SVD 
    % Input 
    % Xs = data matrix 
    % dim1 = number of rows
    % dim2 = number of columns 
    % conv = convergence criterion 
    % iter = maximum number of iterations
    
    % Output 
    % fn = Number of factors/ PCs
    % mse = mean squared error for observed entries compared to the model predictions for these for the optimal number of factors 
    % wmse = winsorized mean squared error for observed entries compared to the model predictions for these for the optimal number of factors
    % R2 = correlation between the observed entries and the model predictions
    % of those entries
    % X_pred = the data matrix as predicted by the optimal number of factors 
    % cutoff = cutoff used to find the number of factors 
    % SVs = singular value matrix found 
    filled_ind = find(~isnan(Xs));
    % find omega 
    if dim1/dim2 ==1
        omega = 2.858; % omega(beta)=2.858 for n*n square matrix
    elseif dim1/dim2 < 1
        beta = dim1/dim2;
        omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
    else
        beta = dim2/dim1; 
        omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
    end 
    [Xfilled, ~] = filldata3(Xs, 'uni',mixtures,conc, 'none', T);
    if dim1<dim2
        fn = dim1;
    else 
        fn=dim2;
    end 
    % find SVs for matrix 
    [U,SVs,D,~,~,~, ~]=missing_svd(Xfilled,fn,1,1,1e-4,2,1);
    % find the threshold 
    y_med = median(diag(SVs)); %ymed = median singular value of the noisy matrix  
    cutoff = omega*y_med; %cutoff= tau= omega(beta)*ymed; matrix
    % Keep modes w/ sig > cutoff; rank chosen as hard cutoff
    fn = length(find(diag(SVs)>cutoff));
    % reconstruct with the new SVs
    SVs(:,fn:end)=0;
     % solve for the correct number of factors 
    [~,~,D,St,X_pred,~,iterations] = missing_svd(Xs,fn,1,1,conv,iter,1, 'avg',mixtures,'none',conc,T); %X,fn,center,scale,conv,max_iter, use_missing,fillmethod,mixtures,whichX,conc,T
    msefill = (sum((X_pred(filled_ind)-Xs(filled_ind)).^2))/length(filled_ind);
    wmsefill = find_wmse_error((Xs(filled_ind)- X_pred(filled_ind)), length(filled_ind));
    RAD = (abs((X_pred(filled_ind)-Xs(filled_ind))./Xs(filled_ind)));
    RAD = RAD(find(~isnan(RAD)));
    RAD = RAD(find(~isinf(RAD)));
    RAD = sum(RAD)/length(RAD);
end 