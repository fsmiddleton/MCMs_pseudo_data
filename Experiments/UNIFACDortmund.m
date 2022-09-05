%% Prediction of excess enthalpy using UNIQUAC parameters  

%FS Middleton 2022/06/20

%% Import data 
% UNIQUAC/ UNIFAC parameters estimated by Aspen using UNIFAC 
clc
clear
filename = 'UNIFACParams.xlsx';

%Table = table of r and q values for each component - structural parameters
Table = readtable(filename,'Sheet', 'RQParams'); 
Table2 = readtable(filename,'Sheet', 'EParams'); 
Table3 = readtable(filename,'Sheet', 'Groups',  'ReadVariableNames',false);
Table4 = readtable(filename,'Sheet', 'RQParamsAll');
fgroup = Table{:,1};
chainlength = Table{:,2};
R = Table{:,3};
Q = Table{:,4};
QOriginal = Table4{:,6};
Eparams = table2array(Table2);
Eparams(isnan(Eparams))=0;

%populate components as numericals 
func_groups = {'Alkane', 'Primaryalcohol', 'Secondaryalcohol','Isoalkanol', 'Tertiaryalcohol','Benzene', 'Toluene', 'Ketone', 'Ketone3','Alkene','Cycloalkane', 'Ester1', 'Ester2','Ester3','Ester4','Ester5','Estercyc', 'Amine', 'Aniline', 'Benzylamine', 'Acid', 'Aldehyde'};
f_num = zeros(length(fgroup),1);
f = 0;
for func = func_groups
    f=f+1;
    indices1 = find(strcmp(func,fgroup));
    f_num(indices1) = f;
end 
components = [f_num chainlength];

fgroup2 = Table3{:,1};
f_num2 = zeros(length(fgroup),1);
f = 0;
for func = func_groups
    f=f+1;
    indices1 = find(strcmp(func,fgroup2));
    f_num2(indices1) = f;
end 
GroupsArray = zeros(size(Table3));
nogroups = size(Table3,1);
GroupsArray(2:end,1) = f_num2(2:end);
GroupsArray(2:end,2) = Table3{2:end,2};
GroupsArray(:,3:end) = Table3{:,3:end};
GroupsArray(isnan(GroupsArray))=0;

sumofgroups1=0;
% find compound 1

% find compound 2
%sumofgroups = 
%% Find Ge/RT for each possible combination of mixtures at every temperature 

%Choose temperatures
Temps =[243.15; 253.15; 263.15; 273.15; 283.15; 288.15; 290.15; 293.15; 296.15; 298.15; 303.15;307.5;309.5; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15; 363.15];
%
conc_interval = 0.1:0.1:0.9;% 0.01:0.01:0.99;

for T = 1:length(Temps)
    T=Temps(T);
    %large concentration interval for ge/RT
    
    mixture = zeros(((size(components,1)^2-size(components,1))/2),4);
    he = zeros(((size(components,1)^2-size(components,1))/2),length(conc_interval));
    
    mixind = 0;
    for component1 = 1:length(components)
        %c1 = index of component in the list of RQParams
        %comp1 = the component functional group and chainlenth, as listed in
        %the data collected 
        comp1 = components(component1,:);
        for component2 =1:length(components)
            comp2 = components(component2,:);
            %check if mixture exists already and if it is a mixture, not a pure
            %component 
            if ~ismember(mixture,[comp1 comp2],'rows') & ~ismember(mixture,[comp2 comp1],'rows') & component1~=component2
                mixind=mixind+1;
                %populate new mixture
                mixture(mixind, :)= [comp1 comp2];
                % find params per compound
                % groups contained in each compound (1 and 2)
                
                r1 = R(component1);
                r2 = R(component2);
                q1 = Q(component1);
                q2 = Q(component2);
                groupsv1 = GroupsArray(component1+1,3:end); 
                groupsv2 = GroupsArray(component2+1,3:end);
                % per group - only for groups, add 1 to the index to ignore row of numbers of groups  
                % returns the indices of the groups required 
                
                countgroup1 = find(GroupsArray(component1+1,3:end)~=0); % groups present in compound 1 
                v1 = GroupsArray(component1+1,countgroup1); % amount of each group
                countgroup2 = (find(GroupsArray(component2+1,3:end)~=0));
                v2 = GroupsArray(component2+1,countgroup2);
                countgroups =[countgroup1 countgroup1];
                groups = {countgroup1 ,countgroup2};
                group1 = GroupsArray(1,countgroup1+2);
                group2 = GroupsArray(1,countgroup2+2);
                %mixture 
                concind = 0;
                for conc = conc_interval
                    concind = concind+1;
                    %find params for each group
                    v1 = GroupsArray(component1+1,3:end); 
                    v2 = GroupsArray(component2+1,3:end);
                    Xmdenom = sum(v1*conc)+sum(v2*(1-conc));
                    Xm1denom = sum(v1*conc);
                    Xm2denom = sum(v2*(1-conc));
                    Xm = (v1*conc+v2*(1-conc))/Xmdenom;
                    Xm1 = (v1*conc)/Xm1denom;
                    Xm2 = v2*(1-conc)/Xm2denom;
                    Qm = (QOriginal(GroupsArray(1,3:end)));
                    groupGC = GroupsArray(1,3:end);
                    XmQm1 = Xm1.*Qm';
                    XmQm2 = Xm2.*Qm';
                    XmQmmix = Xm.*Qm';
                    thetamix = XmQmmix/sum(XmQmmix);
                    theta1 = XmQm1/sum(XmQm1);
                    theta2 = XmQm2/sum(XmQm2);

                    % find a b and c
                    am1 = zeros(length(groupGC)); % mixture 
                    am2 = zeros(length(groupGC)); % mixture 
                    a1 = zeros(length(countgroup1),length(groupGC));% pure substance 1
                    a2 = zeros(length(countgroup2),length(groupGC));% pure substance 2
                    bm1 = zeros(length(groupGC)); % mixture 
                    bm2 = zeros(length(groupGC)); % mixture 
                    b1 = zeros(length(countgroup1),length(groupGC));% pure substance 1
                    b2 = zeros(length(countgroup2),length(groupGC));% pure substance 2
                    cm1 = zeros(length(groupGC)); % mixture 
                    cm2 = zeros(length(groupGC)); % mixture 
                    c1 = zeros(length(countgroup1),length(groupGC));% pure substance 1
                    c2 = zeros(length(countgroup2),length(groupGC));% pure substance 2
                    %find all energetic params, a,b,c, in good orders  
                    for ind = 1:length(countgroup1)
                        for ind2 = 1:length(groupGC)
                            if countgroup1(ind)~=groupGC(ind2) % not the same two groups 
                                index = find(ismember(Eparams(:,1:2), [countgroup1(ind) groupGC(ind2)],'rows'));
                                
                                if index
                                    % row = group of interest; column =
                                    % other groups
                                    am1(ind,ind2) = Eparams(index,3);
                                    bm1(ind,ind2) = Eparams(index,4);
                                    cm1(ind,ind2) = Eparams(index,5);
                                    am1(ind2,ind) = Eparams(index,6);
                                    bm1(ind2,ind) = Eparams(index,7);
                                    cm1(ind2,ind) = Eparams(index,8);
                                    a1(ind,ind2) = Eparams(index,3);
                                    a1(ind+1,ind2) = Eparams(index,6);                              
                                    b1(ind,ind2) = Eparams(index,4);
                                    b1(ind+1,ind2) = Eparams(index,7);
                                    c1(ind,ind2) = Eparams(index,5);
                                    c1(ind+1,ind2) = Eparams(index,8);
                                end 
                            end 
                        end
                    end
                    for ind = 1:length(countgroup2)
                        for ind2 = 1:length(groupGC)
                            if countgroup2(ind)~=groupGC(ind2) % not the same two groups 
                                index = find(ismember(Eparams(:,1:2), [countgroup2(ind) groupGC(ind2)],'rows'));
                                forwards = 1;
                                
                                if index
                                    am2(ind,ind2) = Eparams(index,3);
                                    bm2(ind,ind2) = Eparams(index,4);
                                    cm2(ind,ind2) = Eparams(index,5);
                                    am2(ind2,ind) = Eparams(index,6);
                                    bm2(ind2,ind) = Eparams(index,7);
                                    cm2(ind2,ind) = Eparams(index,8);
                                    a2(ind,ind2) = Eparams(index,3);
                                    a2(ind+1,ind2) = Eparams(index,6);
                                    b2(ind,ind2) = Eparams(index,4);
                                    b2(ind+1,ind2) = Eparams(index,7);
                                    c2(ind,ind2) = Eparams(index,5);
                                    c2(ind+1,ind2) = Eparams(index,8);
                                end 
                            end 
                        end 
                    end 
                    % phi calcs 
                    psimix1 = exp(-1/T*(am1+bm1*T+cm1*T^2));
                    psimix2 = exp(-1/T*(am2+bm2*T+cm2*T^2));
                    psi1 = exp(-1/T*(a1+b1*T+c1*T^2));
                    psi2 = exp(-1/T*(a2+b2*T+c2*T^2));
                    %have all the params, calculate values 
                    residual = zeros(2);
                    for i = 1:2 % compounds in the mixture 
                        residual(i)=0;
                        
                        dTm = zeros(length(groups{i}),2);
                        dTs = zeros(length(groups{i}),2);
                        for k = 1:length(groups{i}) % groups in the compound
                            %mixture 
                            
                            
                            if i==1
                                psim = psimix1(:,k)';
                                psik = psimix1(k,:);
%                                 thetamix = thetamix(:,k)';
                         
                                bkm = bm1(:,k)';
                                bkms = b1(k,:);
                                ckm = cm1(:,k)';
                                ckms = c1(k,:);
                                %mixture
                                dTm(k,i)=Qm(k)*sum((thetamix.*((bkm+log(psim)+2*ckm*T).*psim)/sum(thetamix.*psik)) + (thetamix*sum(psim.*thetamix)./((sum(psim.*thetamix))^2).*(bkm-b1(k,:)+log(psik./psim)+(ckm-c1(k,:))*2*T)));
                                %pure substance 
                                dTs(k,i)=Qm(k)*sum(theta1.*((b1(k,:)+log(psi1(k,:))+2*c1(k,:)*T).*psi1(k,:))/sum(theta1.*psi1(k,:)) + theta1*sum(psi1(k,:).*theta1)./((sum(psi1(k,:).*theta1))^2).*(bkms-b1(k,:)+log(psi1(k,:)./psim)+(ckms-c1(k,:))*2*T));
                            else
                                psim = psimix2(:,k)';
                                psik = psimix2(k,:);
                                bkm = bm2(:,k)';
                                bkms = b2(k,:);
                                ckm = cm2(:,k)';
                                ckms = c2(k,:);
                                dTm(k,i)=Qm(k)*sum(thetamix.*((bkm+log(psim)+2*ckm*T).*psim)/sum(thetamix.*psim) + (thetamix.*sum(psim.*thetamix)./((sum(psim.*thetamix))^2).*(bkm-b2(k,:)+log(psik./psim)+(ckm-c2(k,:))*2*T)));
                                %pure substance 
                                dTs(k,i)=Qm(k)*sum(theta2.*((b2(k,:)+log(psi2(k,:))+2*c2(k,:)*T).*psi2(k,:))./sum(theta2.*psi2(k,:)) + theta2.*sum(psi2(k,:).*theta2)./((sum(psi2(k,:).*theta2))^2).*(bkms-b2(k,:)+log(psi2(k,:)./psim)+(ckms-c2(k,:))*2*T));
                            end
                        end 
                        residual(i) = sum(groups{i}.*(dTm(:,i)-dTs(:,i))'); % groups is a row vector 
                    end
                    % calculate he and finally output it 
                    he(mixind,concind) = -8.314*T*(conc*residual(1)+(1-conc)*residual(2));
                end
            end 
        end 
    end 
     filenametemp = strcat('heUNIFACforT=', num2str(T), '.mat');
     save(filenametemp);
end 