%% Find scree plots and hard-thresholded ranks (Gavish)
% FS Middleton 12 November 2022

%% Import .m file with 4-way array
clc
clear
% decide on which interval and how many temperatures to evaluate
interval = 0.05;
NoT = 3;
load(strcat('HE4wayArrayPolyAll',num2str(NoT),'.mat'))
conc_interval = interval:interval:(1-interval);
X = HE_data_sparse(:,:,1:length(conc_interval),:);
Xsign = sign(X);
Xscale = Xsign.*log(Xsign.*X);
dim = size(X);
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end 
[compoundnames,codesout] = findnames(comps);
%% Rank determined from SVs
concentrations = conc_interval;
fn = min([dim(2),dim(4)]);
y1 = zeros(fn,length(conc_interval));
y2 = zeros(fn,length(conce_interval));

%necessary for plots
indexplot=fn;
whichX = 'none';
if strcmp(whichX,'scale')
    Xss=Xscale;
elseif  strcmp(whichX,'sign')
    Xss=Xsign;
else
    Xss=X;
end
i=0;

for c = (conc_interval(1:9))
    i=i+1;
    f = figure(i);
    clf
    title(strcat('Concentration ',num2str(conc_interval(i))))
    for m = 1:dim(1)
        Xs=reshape(Xss(m,:,i,:),dim2,dim(4));
        %SVD
        Xf = filldata3(Xs, 'uni',mixtures,c,whichX,T);
        [~,D,~,~,~, ~]=missing_svd(Xf,fn,1,1,1e-8,1,0);
        %Plot
        y1(:,i) = diag(D);
        y2(:,i) = cumsum(y1(:,i))/sum(y1(:,i));
        subplot(2,1,1)
        xlabel('Number of factors')
        ylabel('Singular Values')
        plot(1:indexplot, y1(1:indexplot,i))
        hold on
        subplot(2,1,2)
        xlabel('Number of factors')
        ylabel('Cumulative contribution')
        plot(1:indexplot,y2(1:indexplot,i))
        hold on
    end 
    title(strcat('Concentration ',num2str(conc_interval(i))))
    legend((compoundnames), 'Location', 'southeast','FontSize',6,'TextColor','black','NumColumns',4)
end 
hold off 

%% Gavish
clc
clear

NoT =7;
interval = 0.05;
filename = strcat('HE4wayArrayPolyAll',num2str(NoT),'.mat');
load(filename, 'HE_data_sparse',  'comps', 'mixtureT','mixture', 'Temps')
conc_interval = interval:interval:(1-interval);
X = HE_data_sparse(:,:,1:length(conc_interval),:);
Xsign = sign(X);
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end 
interval = 0.05;
conc_interval = interval:interval:(1-interval);
indend = size(HE_data_sparse,1);
X = HE_data_sparse(1:indend,1:indend,1:length(conc_interval),:);
Xsign = sign(X);
Xscale = Xsign.*log(Xsign.*X);
dim = size(X);
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end
whichX = 'none';
if strcmp(whichX,'scale')
    Xss=Xscale;
elseif  strcmp(whichX,'sign')
    Xss=Xsign;
else
    Xss=X;
end

j=0;

Xs = filldata3(Xss,'uni',mixtures,conc_interval,whichX,Temps);

for c = conc_interval % for concentrations
    j=j+1; 
    for m = 1:dim1 %for compounds
        Xfilled = reshape(Xs(m,:,j,:),dim2,[]);
        [min_fn(j,m),minmse_fill(j,m),minwmse_fill(j,m),minRAD(j,m), X_pred_best,~,~,iter(j,m)] = solveGavish(Xfilled, dim2, dim(4),1e-10,2, mixtures, c,Temps);
    end 
end

