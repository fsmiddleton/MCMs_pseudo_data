function [compoundnames,codesout] = findnames(codesin)
    % Find the names of the compounds in mixtures, given the codes of the
    % compounds
    % Input
    % codesin = array of codes with two columns indicating functional group
    % and chain length 
    % 
    % Outputs 
    % compoundnames = names of the compounds 
    % codesout = codes in the same order as the names outputted 
    compounds = {'Methane','Ethane', 'Propane', 'Butane', 'Pentane', 'Hexane', 'Heptane', 'Octane', 'Nonane', 'Decane', 'Dodecane', 'Methanol', 'Ethanol', '1-Propanol', '1-Butanol', '1-Pentanol', '1-Hexanol', '1-Heptanol', '1-Octanol', '1-Nonanol', '1-Decanol', '2-Propanol', '2-Butanol', '2-Pentanol', '2-Hexanol','2-Octanol', 'Isobutanol', 'Tertbutylalcohol', 'Cyclopentane', 'Cyclohexane','Cycloheptane', 'Cyclooctane','Benzene', 'Toluene', 'Ethanal', 'Propanal', 'Butanal', 'Pentanal', 'Hexanal', 'Ethene', '1-Pentene', '1-Hexene', '1-Heptene', '1-Octene', '1-Nonene', '1-Decene', 'Acetone', 'Butanone', '2-Propanone', '2-Hexanone', '2-Heptanone', '2-Octanone','Formic acid', 'Acetic acid', 'Propionic acid', 'Butyric acid', 'Pentanoic acid', 'Hexanoic acid', 'Heptanoic acid', 'Octanoic acid', 'Nonanoic acid', 'Decanoic acid', '3-Pentanone', '3-Hexanone', '3-Heptanone', '3-Octanone', 'Methyl formate', 'Methyl acetate', 'Methyl propionate', 'Methyl butyrate', 'Methyl pentanoate', 'Methyl hexanoate', 'Methyl benzoate', 'Ethyl benzoate', 'Ethyl formate', 'Ethyl acetate', 'Ethyl propionate', 'Ethyl butyrate', 'Ethyl pentanoate', 'Ethyl hexanoate', 'Ethyl heptanoate', 'Ethyl octanoate' , 'Propyl formate','Propyl acetate', 'Propyl propionate', 'Propyl butyrate', 'Butyl formate', 'Butyl acetate', 'Butyl butyrate', 'Pentyl acetate', 'Propylamine', 'Butylamine', 'Pentyalamine', 'Hexylamine', 'Heptylamine', 'Octylamine', 'Nonylamine', 'Decylamine', 'Aniline', 'Benzylamine'};
    codes = [ones(11,1), ([1:10,12])'; 2*ones(10,1), (1:10)'; 3*ones(5,1), ([3:6,8])'; 4, 4; 5, 4;11*ones(4,1), (5:8)'; 6, 0; 7, 0; 22*ones(5,1), (2:6)'; 10*ones(7,1), ([2,5:10])'; 8*ones(6,1), (3:8)'; 21*ones(10,1), (1:10)'; 9*ones(4,1), (5:8)';12*ones(6,1), (1:6)'; 17*ones(2,1), (1:2)'; 13*ones(8,1), (1:8)'; 14*ones(4,1), (1:4)'; 15*ones(3,1), ([1,2,4])'; 16, 1; 18*ones(8,1), (3:10)'; 19, 0; 20, 0];
    indices = find(ismember(codes, codesin,'rows'));
    for i = 1:length(indices)
        compoundnames{i} = compounds{indices(i)};
        codesout(i,:) =codes(indices(i),:);
    end 
end 
 