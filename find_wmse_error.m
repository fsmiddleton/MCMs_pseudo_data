function [wmse]=find_wmse_error(errors, count)
    errors = reshape(errors,prod(size(errors)),1);
    perc5=prctile(errors,5,'all');
    perc95=prctile(errors,95,'all');
    errors(errors<perc5)=perc5;
    errors(errors>perc95)=perc95;
    wmse = (sum((errors).^2))/count;
end 
