%% Plotting distributions to find differences between them 
% Francesca Middleton - 2022/03/11

clc
clear 
clf 
% Cauchy functions from Peder Axensten (2022). cauchy (https://www.mathworks.com/matlabcentral/fileexchange/11749-cauchy), MATLAB Central File Exchange. Retrieved March 11, 2022.
%% generate distributions 

x = -3:.1:3;
gauss = normpdf(x,0,1);
cauchy = cauchypdf(x,0,1);

plot(x,gauss)
hold on 
plot(x, cauchy)
legend('Gaussian', 'Cauchy')