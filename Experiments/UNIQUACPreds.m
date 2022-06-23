%% Prediction of excess enthalpy using UNIQUAC parameters  

%FS Middleton 2022/06/20

%% Import data 

filename = 'UNIQUACparams.xlsx';

%Table = table of r and q values for each component 
Table = readtable(filename,'Sheet', 'RQParams'); 
fgroup = Table{:,1};
chainlength = Table{:,2};
R = Table{:,3};
Q = Table{:,4};
Q2a = Table{:,5};
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

Table2 = readtable(filename,'Sheet', 'UNIQUAC Params'); 
%% Find Ge/RT for each possible combination of mixtures at a chosen temperature 

%Choose temperature 
T = 298.15; %K

%large concentration interval for ge/RT
conc_interval = 0.01:0.01:0.99;
mixture = zeros(4,(length(fgroup)^2/2-length(group)));
mixind = 0;
for c1 = 1:length(components)
    comp1 = components(c1);
    for c2 = 1:length(components)
        comp2 = components(c2);
        %check if mixture exists already and if it is truly a mixture
        if ~ismember(mixture,[comp1 comp2],'rows') && ~ismember(mixture,[comp2 comp1],'rows') && c1~=c2
            mixind=mixind+1;
            %populate new mixture
            mixture(:, mixind)= [comp1 comp2];
            r1 = R(c1);
            r2 = R(c2);
            q1 = Q(c1);
            q2 = Q(c2);
            q1a = Qa(c1);
            q2a = Qa(c2);
            %find a12 and a21
            
            for c = 1:length(conc_interval)
                conc = conc_interval(c);
                phi1 = conc*r1/(conc*r1+(1-conc)*r2);
                phi2 = (1-conc)*r2/(conc*r1+(1-conc)*r2);
                theta1 = conc*q1/(conc*q1+(1-conc)*q2);
                theta2 = (1-conc)*q2/(conc*q1+(1-conc)*q2);
                thetaa1= conc*q1a/(conc*q1a+(1-conc)*q2a);
                thetaa2 = (1-conc)*q2a/(conc*q1a+(1-conc)*q2a);
                tau12
                tau21
                geRTr = -conc*q1a*log(thetaa1+thetaa2*tau21)-(1-conc)*q2a*log(thetaa2+thetaa1*tau12);
                geRTc = conc*log(theta1/conc)+(1-conc)*log(theta2/(1-conc))+5*(conc*q1*log(theta1/phi1)+(1-conc)*q2*log(theta2/phi2));
                geRT = geRTr + geRTc;
            end 
        end 
    end 
end 

% Find he = -T^2(dGe/T/dT) for this temperature 