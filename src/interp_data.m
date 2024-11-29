function [HE, uncertainty, orderPolyfit, HEpredout, errorpredout, R2]=interp_data(data, conc_interval)
    % Interpolate data for use in the matrices
    %
    % Inputs 
    % data = the whole table of data to be interpolated
    % conc_interval = list of concentrations for interpolation
    %
    % Outputs  
    % HE = excess enthalpy data at each concentration value for the mixture
    % uncertainty = uncertainty associated with each HE value 
    % orderPolyfit = order of the fit for each mixture 
    % HEpredout = predictions of HE at the same concentrations as the
    % imported data
    % erorpredout = errors of the HEpredout predictions
    % R2 = fit of the data 
    
    HE_original = data.Excessenthalpy;
    comp = data.Compositioncomponent1;
    % remove the 0 and 1 data points 
    ind_keep = find(HE_original ~=0);
    HE_new= HE_original(ind_keep);
    comp_new = comp(ind_keep);
    numberofpoints = length(comp_new);
    if ind_keep
        comp_new(numberofpoints+1)=0;
        comp_new(numberofpoints+2)=1;
        HE_new(numberofpoints+1)=0;
        HE_new(numberofpoints+2)=0;
    end 
    if length(comp_new)==3 
        % only one entry at 50 mole % - the rest of the data is not filled
        if comp_new(1) == 0.5
            HE = nan(length(conc_interval),1);
            HE(find(conc_interval==0.5))=HE_new(1);
            orderPolyfit = 1;
            uncertainty =zeros(length(HE),1);
            HEpredout = zeros(300,1);
            HEpredout(1:length(HE),1)=HE;
            errorpredout = zeros(300,1);
            R2 = 1;
        else 
            % only 1 entry not at 50%
            % round to the nearest concentration by fitting a second order
            % polynomial and using the value generated at the nearest
            % concentration 
            concdiff = conc_interval-comp_new(1);
            index = find(concdiff==min(abs(concdiff)));
            [p,S] = polyfit(comp_new, HE_new,2);
            [HEtemp,uncertainty] = polyval(p,conc_interval, S);
            HE = nan(length(HEtemp),1);
            HE(index) = HEtemp(index);
            [HEpred, errorpred] = polyval(p, comp_new,S);
            HEpredout = zeros(300,1);
            HEpredout(1:length(HEpred),1)=HEpred;
            errorpredout = zeros(300,1);
            errorpredout(1:length(errorpred),1)=errorpred;
            R2 = 1;
            orderPolyfit = 1;
        end 
    else 
        % check for data that does not have many points outside a certain interval 
        maxdiff = max(comp)-min(comp);
        %interpolate the data 
        if maxdiff <0.3 || numberofpoints <3
            %just fit a parabola to the small amount of data 
            maxOrder = 2;
        else 
            maxOrder = 7;
        end 
        error = zeros(maxOrder-1,1);
        ind = 1;
        for i =2:maxOrder
            [p,S,mu] = polyfit(comp_new, HE_new,i);
            [~,uncertainty] = polyval(p,conc_interval, S,mu);
            error(ind) = sum(uncertainty)/length(uncertainty);
            ind = ind+1;
        end 
        %populate the data 
        orderPolyfit = find(error == min(error))+1;
        [p,S] = polyfit(comp_new, HE_new,orderPolyfit);
        [HE,uncertainty] = polyval(p,conc_interval, S);
        [HEpred, errorpred] = polyval(p, comp_new,S);
        HEpredout = zeros(300,1);
        HEpredout(1:length(HEpred),1)=HEpred;
        errorpredout = zeros(300,1);
        errorpredout(1:length(errorpred),1)=errorpred;
        R2 = corrcoef(HEpred, HE_new);
        R2 = R2(2,1);
    end 
end 