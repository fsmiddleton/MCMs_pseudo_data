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
Qgroup = Table4{:,6};
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
%% Find He for each possible combination of mixtures at every temperature 

%Choose temperatures
Temps =298.15; %[243.15; 253.15; 263.15; 273.15; 283.15; 288.15; 290.15; 293.15; 296.15; 298.15; 303.15;307.5;309.5; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15; 363.15];
%
conc_interval = .5; %0.05:0.05:0.95;% 0.01:0.01:0.99;

for T = (Temps)
    
    mixture = zeros(((size(components,1)^2-size(components,1))/2),4);
    he = zeros(((size(components,1)^2-size(components,1))/2),length(conc_interval));
    
    mixind = 0;
    for component1 = 6 %1:length(components)
        %c1 = index of component in the list of RQParams
        %comp1 = the component functional group and chainlenth, as listed in
        %the data collected 
        comp1 = components(component1,:);
        for component2 =17 %1:length(components)
            comp2 = components(component2,:);
            %check if mixture exists already and if it is a mixture, not a pure
            %component 
            if ~ismember(mixture,[comp1 comp2],'rows') & ~ismember(mixture,[comp2 comp1],'rows') & component1~=component2
                mixind=mixind+1;
                %populate new mixture
                mixture(mixind, :)= [comp1 comp2];
                % find params per compound
                % groups contained in each compound (1 and 2)
                
%                 r1 = R(component1);
%                 r2 = R(component2);
%                 q1 = Q(component1);
%                 q2 = Q(component2);
                
                % Extract parameters 
                 % each component's row vector of vs
                groupsv1 = GroupsArray(component1+1,3:end);
                groupsv2 = GroupsArray(component2+1,3:end);
                % indices and the number of each group
                indgroup1 = find(GroupsArray(component1+1,3:end)~=0); % groups present in compound 1
                groupNo1 = GroupsArray(1, 2+indgroup1);
                indgroup2 = (find(GroupsArray(component2+1,3:end)~=0));
                groupNo2 = GroupsArray(1, 2+indgroup2);
                %both sets of groups 
                groupsNo = union(groupNo1, groupNo2);
                indgroups = union(indgroup1,indgroup2);
                % pure compound params
                v1 = GroupsArray(component1+1,2+indgroup1); % number of each group
                v2 = GroupsArray(component2+1,2+indgroup2);
                Q1 = Qgroup(groupNo1);
                Q2 = Qgroup(groupNo2);
                %group params 
                Q12 = Qgroup(groupsNo);
                v1m = GroupsArray(component1+1,2+indgroups); % amount of each group
                v2m = GroupsArray(component2+1,2+indgroups);
                %interaction parameters 
                % individual compounds 
                % define the parameters, then fill in a for loops
                a1 = zeros(length(groupNo1));
                b1 = zeros(length(groupNo1));
                c1 = zeros(length(groupNo1));
                a2 = zeros(length(groupNo2));
                b2 = zeros(length(groupNo2));
                c2 = zeros(length(groupNo2));
                am = zeros(length(groupsNo));
                bm = zeros(length(groupsNo));
                cm = zeros(length(groupsNo));
                % compound 1 
                for k = 1:length(groupNo1)
                    for m = 2:length(groupNo1)
                        group1 = groupNo1(k);
                        group2 = groupNo1(m);
                        index = find(ismember(Eparams(:,1:2), [group1 group2],'rows'));
                        if index
                            a1(k,m) = Eparams(index,3);
                            a1(m,k) = Eparams(index,6);
                            b1(k,m) = Eparams(index,4);
                            b1(m,k) = Eparams(index,7);
                            c1(k,m) = Eparams(index,5);
                            c1(m,k) = Eparams(index,8);
                        end 
                    end 
                end 
                %compound 2
                for k = 1:length(groupNo2)
                    for m = 2:length(groupNo2)
                        group1 = groupNo2(k);
                        group2 = groupNo2(m);
                        index = find(ismember(Eparams(:,1:2), [group1 group2],'rows'));
                        if index
                            a2(k,m) = Eparams(index,3);
                            a2(m,k) = Eparams(index,6);
                            b2(k,m) = Eparams(index,4);
                            b2(m,k) = Eparams(index,7);
                            c2(k,m) = Eparams(index,5);
                            c2(m,k) = Eparams(index,8);
                        end
                    end 
                end 
                %mixture 
                for k = 1:length(groupsNo)
                    for m = 2:length(groupsNo)
                        group1 = groupsNo(k);
                        group2 = groupsNo(m);
                        index = find(ismember(Eparams(:,1:2), [group1 group2],'rows'));
                        if index
                            am(k,m) = Eparams(index,3);
                            am(m,k) = Eparams(index,6);
                            bm(k,m) = Eparams(index,4);
                            bm(m,k) = Eparams(index,7);
                            cm(k,m) = Eparams(index,5);
                            cm(m,k) = Eparams(index,8);
                        end 
                    end 
                end 
                psi1 = exp(-1/T*(a1+b1.*T+c1*T^2));
                psi2 = exp(-1/T*(a2+b2.*T+c2*T^2));
                psim = exp(-1/T*(am+bm.*T+cm*T^2));
                
                %mixture 
                concind = 0;
                for conc = conc_interval
                    concind = concind+1;
                    % Find Xm and thetam for this concentration 
                    Xmdenom = sum(v1m*conc)+sum(v2m*(1-conc));
                    Xm1denom = sum(v1*conc);
                    Xm2denom = sum(v2*(1-conc));
                    Xm = [conc 1-conc]*[v1m;v2m]/Xmdenom;
                    Xm1 = (v1*conc)/Xm1denom;
                    Xm2 = v2*(1-conc)/Xm2denom;
                    XmQm1 = Xm1.*Q1';
                    XmQm2 = Xm2.*Q2';
                    XmQmmix = Xm.*Q12';
                    thetam = XmQmmix/sum(XmQmmix);
                    theta1 = XmQm1/sum(XmQm1);
                    theta2 = XmQm2/sum(XmQm2);


                    %have all the params, calculate values 
                    residual = zeros(2);
                    % compound1
                    
                    dTp = zeros(length(v1),1);
                    dTm = zeros(length(v1),1);
                    % pure compound 
                    for k=1:length(v1)%k
                        %for each group in compound 1, calculate the terms
                        %pure compound 
                        indices = 1:length(v1);
                        mindex = indices(indices~=k);
                        y = zeros(length(mindex),1);
                        if mindex
                            for m=1:length(mindex)%m
                                nindex = indices(indices~=m);%n
                                % pure compound 
                                y(m)= [psi1(m,k) psi1(k,m)]*[(b1(m,k)+log(psi1(m,k))+2*c1(m,k)*T)/sum(theta1(nindex).*psi1(nindex,m)') ; psi1(k,m)/sum((theta1(nindex).*psi1(nindex,m)').^2)*sum(theta1(nindex).*psi1(nindex,m)'.*(b1(k,m)-b1(nindex,m)'+log(psi1(k,m)./psi1(nindex,m)')+2*T.*(c1(k,m)-c1(nindex,m)')))]; 
                            end 
                            dTp(k) = v1(k)*Q1(k)*theta1(mindex)*y;
                        end 
                        %mixture            
                        indices = 1:length(v1m);
                        mindex = indices(indices~=k);
                        %mixture 
                        y = zeros(length(mindex),1);
                        for m=1:length(mindex)%m
                            nindex = indices(indices~=m);%n
                            y(m)= [psim(m,k) psim(k,m)]*[(bm(m,k)+log(psim(m,k))+2*cm(m,k)*T)/sum(thetam(nindex).*psim(nindex,m)') ; psim(k,m)/sum((thetam(nindex).*psim(nindex,m)').^2)*sum(thetam(nindex).*psim(nindex,m)'.*(bm(k,m)-bm(nindex,m)'+log(psim(k,m)./psim(nindex,m)')+2*T.*(cm(k,m)-cm(nindex,m)')))]; 
                        end
                        dTm(k) = v1(k)*Q1(k)*thetam(mindex)*y;
                    end 
                    %calculate term for compound 1
                    residual1= dTm-dTp;
                    
                    %compound 2
                    dTp = zeros(length(v2),1);
                    dTm = zeros(length(v2),1);
                    % pure compound 
                    for k=1:length(v2)%k
                        %for each group in compound 1, calculate the terms
                        indices = 1:length(v2);
                        mindex = indices(indices~=k);
                        y = zeros(length(mindex),1);
                        if mindex
                            for m=1:length(mindex)%m
                                nindex = indices(indices~=m);%n
                                y(m)= [psi2(m,k) psi2(k,m)]*[(b2(m,k)+log(psi2(m,k))+2*c2(m,k)*T)/sum(theta2(nindex).*psi2(nindex,m)') ; psi2(k,m)/sum((theta2(nindex).*psi2(nindex,m)').^2)*sum(theta2(nindex).*psi2(nindex,m)'.*(b2(k,m)-b2(nindex,m)'+log(psi2(k,m)./psi2(nindex,m)')+2*T.*(c2(k,m)-c2(nindex,m)')))]; 
                            end 
                            dTp(k) = v2(k)*Q2(k)*theta2(mindex)*y;
                        end 
                        
                        %mixture 
                        indices = 1:length(v2m);
                        mindex = indices(indices~=k);
                        %mixture 
                        y = zeros(length(mindex),1);
                        for m=1:length(mindex)%m
                            nindex = indices(indices~=m);%n
                            y(m)= [psim(m,k) psim(k,m)]*[(bm(m,k)+log(psim(m,k))+2*cm(m,k)*T)/sum(thetam(nindex).*psim(nindex,m)') ; psim(k,m)/sum((thetam(nindex).*psim(nindex,m)').^2)*sum(thetam(nindex).*psim(nindex,m)'.*(bm(k,m)-bm(nindex,m)'+log(psim(k,m)./psim(nindex,m)')+2*T.*(cm(k,m)-cm(nindex,m)')))]; 
                        end
                        dTm(k) = v2(k)*Q2(k)*thetam(mindex)*y;
                    end 
                    residual2 = dTm-dTp;
                    % calculate he and output to correct mixture and
                    % concentration index 
                    he(mixind,concind) = -8.314*T*(conc*sum(residual1)+(1-conc)*sum(residual2));
                end
            end 
        end 
    end 
     filenametemp = strcat('heUNIFACforT=', num2str(T), '.mat');
     %save(filenametemp);
end 