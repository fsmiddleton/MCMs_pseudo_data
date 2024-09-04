function [wmse]=find_wmse_error(errors, count)
    % Find the wMSE of the errors inputted to this function
    % Inputs 
    % errors = array of errors
    % count = count of the number of entries used in errors
    %
    % Output 
    % wmse = wMSE of the errors
    
    errors = reshape(errors,numel(errors),1);
    perc5=prctile(errors,5,'all');
    perc95=prctile(errors,95,'all');
    errors(errors<perc5)=perc5;
    errors(errors>perc95)=perc95;
    wmse = (sum((errors).^2))/count;
end 
