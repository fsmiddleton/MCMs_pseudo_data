clear 
clc 
clf

xC = [2;.3;];
sig = [2;1;];
theta=pi/3;
R=[cos(theta) -sin(theta);
    sin(theta) cos(theta)];
nPoints = 10000;
X = R*diag(sig)*randn(2,nPoints)+diag(xC)*ones(2,nPoints);

subplot(1,2,1)
scatter(X(1,:),X(2,:),'k','LineWidth',0.5)
hold on, box on, grid on
axis([-6 8 -6 8])

Xavg= mean(X,2);
B=X-Xavg*ones(1,nPoints);% mean subtracted data 
[U,S,V]=svd(B/sqrt(nPoints),'econ');%PCAs

subplot(1,2,2)
scatter(X(1,:),X(2,:),'k','LineWidth',0.5)
hold on, box on, grid on
axis([-6 8 -6 8])
theta =(0:0.01:1)*2*pi;
%1std confidence interval 
Xstd=U*S*[cos(theta); sin(theta)];
%predictions with 1, 2, and 3 std devs
plot(Xavg(1)+Xstd(1,:), Xavg(2)+Xstd(2,:), 'r', 'Linewidth',1.5)
plot(Xavg(1)+2*Xstd(1,:), Xavg(2)+2*Xstd(2,:), 'r', 'Linewidth',1.5)
plot(Xavg(1)+3*Xstd(1,:), Xavg(2)+3*Xstd(2,:), 'r', 'Linewidth',1.5)

%% Ovarian cancer 
clf
clc
clear
close all
%%
load ovariancancer;
data=obs;
[U,S,V]=svd(data,'econ');

figure
subplot(1,2,1)
semilogy(diag(S), 'k-o','LineWidth',1)
set(gca,'FontSize',15), axis tight, grid on 
subplot(1,2,2)
plot(cumsum(diag(S))/sum(diag(S)),'k-o','LineWidth',1)
set(gca,'FontSize',15), axis tight, grid on 
%set(gcf, 'Position', [1400,100,3*600,3*250])
%%

figure, hold on
for i=1:size(data,1)
    x=V(:,1)'*data(i,:)';
    y=V(:,2)'*data(i,:)';
    z=V(:,3)'*data(i,:)';
    if (grp{i}=='Cancer')
        plot3(x,y,z,'rx','LineWidth',3);
    else 
        plot3(x,y,z, 'bo','LineWidth',2);
    end 
end 
xlabel('PC1'), ylabel('PC2'), zlabel('PC3')
view(85,25), grid on , set(gca,'FontSize',15)
%set(gcf, 'Position', [1400,100,1200,900])

