clear 
clc 
clf

xC = [2;0.5;];
sig = [2;.5;];
theta=pi/3;
R=[cos(theta) -sin(theta);
    sin(theta) cos(theta)];
nPoints = 1000;
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


hold off
