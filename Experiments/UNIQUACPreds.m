%% Prediction of excess enthalpy using UNIQUAC parameters  

%FS Middleton 2022/06/20

%% Import data 
clc
clear
filename = 'UNIQUACparams.xlsx';

%Table = table of r and q values for each component 
Table = readtable(filename,'Sheet', 'RQParams'); 
fgroup = Table{:,1};
chainlength = Table{:,2};
R = Table{:,3};
Q = Table{:,4};
Qa = Table{:,5};
AspenName = Table{:,6};
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

Tablereal = readtable(filename,'Sheet', 'UNIQUACreal', 'EmptyValue', 0); 
Tableest = readtable(filename,'Sheet', 'UNIQUACest', 'EmptyValue', 0); 
Tablenames = readtable(filename, 'Sheet', 'AspenNames');
AspenNames = Tablenames{:,1};
RealNames1 = Tablereal{:, 1};
RealNames2 = Tablereal{:,2};
RealParams = Tablereal{:,5:8};
EstNames = strcat(Tableest{:,1}, Tableest{:,2});
EstParams = Tableest{:,3:6};
AspenComps = Tablenames{:,[3,4]};
%% Find Ge/RT for each possible combination of mixtures at every temperature 

%Choose temperatures
Temps = [243.15; 253.15; 263.15; 273.15; 283.15; 288.15; 290.15; 293.15; 296.15; 298.15; 303.15; 308.15; 313.15; 318.15; 323.15; 328.15; 333.15; 343.15; 348.15; 353.15; 363.15];
conc_interval = 0.01:0.01:0.99;
for T =1:length(Temps)

    %large concentration interval for ge/RT
    
    mixture = zeros((length(fgroup)^2/2-length(fgroup)),4);
    geRT = zeros((length(fgroup)^2/2-length(fgroup)),length(conc_interval));
    mixind = 0;
    for c1 = 1:length(components)
        %c1 = index of component in the list of RQParams
        %comp1 = the component functional group and chainlenth, as listed in
        %the data collected 
        comp1 = components(c1,:);
        for c2 = 1:length(components)
            comp2 = components(c2,:);
            %check if mixture exists already and if it is a mixture, not a pure
            %component 
            if ~ismember(mixture,[comp1 comp2],'rows') & ~ismember(mixture,[comp2 comp1],'rows') & c1~=c2
                mixind=mixind+1;
                %populate new mixture
                mixture(mixind, :)= [comp1 comp2];
                r1 = R(c1);
                r2 = R(c2);
                q1 = Q(c1);
                q2 = Q(c2);
                q1a = Qa(c1);
                q2a = Qa(c2);
                %find a12 and a21
                name1 = AspenNames(find(ismember(AspenComps,comp1,'rows')));
                name2 = AspenNames(find(ismember(AspenComps,comp2,'rows')));
                namecomp1 = (find(strcmp(name1,RealNames1)));
                namecomp2 = (find(strcmp(name2,RealNames2)));
                ind1 = find(ismember(namecomp1,namecomp2));
                if ~isempty(ind1)
                    %mixture's real UNIQUAC parameters exist in Aspen
                    %literature
                    index = namecomp1(ind1);
                    a12 = RealParams(ind1,1);
                    a21 = RealParams(ind1,2);
                    b12 = RealParams(ind1,3);
                    b21 = RealParams(ind1,4);
                else
                    %Estimated using UNIFAC 
                    namecat = strcat(name1,name2);
                    %does not work 
                    ind = find(strcmp(namecat,EstNames)); % unifac and unifac-dmd exist
                    if isempty(ind)
                        disp('No parameters for this mixture')
                        disp(namecat)
                        ind=1;
                    elseif length(ind)>1
                        ind = ind(2); %Chooses UNIFAC-Dortmund over UNIFAC
                    end 
                    a12 = EstParams(ind,1);
                    a21 = EstParams(ind,2);
                    b12 = EstParams(ind,3);
                    b21 = EstParams(ind,4);
                end

                %calculate taus for this mixture 
                tau12 = exp(-a12/T-b12);
                tau21 = exp(-a21/T-b21);

                for c = 1:length(conc_interval)
                    conc = conc_interval(c);
                    phi1 = conc*r1/(conc*r1+(1-conc)*r2);
                    phi2 = (1-conc)*r2/(conc*r1+(1-conc)*r2);
                    theta1 = conc*q1/(conc*q1+(1-conc)*q2);
                    theta2 = (1-conc)*q2/(conc*q1+(1-conc)*q2);
                    thetaa1= conc*q1a/(conc*q1a+(1-conc)*q2a);
                    thetaa2 = (1-conc)*q2a/(conc*q1a+(1-conc)*q2a);
                    geRTr = -conc*q1a*log(thetaa1+thetaa2*tau21)-(1-conc)*q2a*log(thetaa2+thetaa1*tau12);
                    geRTc = conc*log(theta1/conc)+(1-conc)*log(theta2/(1-conc))+5*(conc*q1*log(theta1/phi1)+(1-conc)*q2*log(theta2/phi2));
                    geRT(mixind, c) = geRTr + geRTc;
                end 
            end 
        end 
    end 
    filenametemp = strcat('geRTforT=', num2str(T), '.mat');
    save(filenametemp);
end 
%% Find he = -T^2(dGe/T/dT) for this temperature 