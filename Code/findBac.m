function [bacnames,bacgroups] = findBac(mixtures)
    % Find the name of the BAC type to which each compound belongs, and the
    % number of the BAC group that a mixture belongs to 
    % Input 
    % mixtures = codes of mixtures, which should have four columns 
    %
    % Outputs 
    % bacnames = BAC type to which each compound belongs
    % bacgroup = BAC group to which each mixture belongs 
    bacNamegroups = {'NA', 'HA', 'SA', 'HD'};
    
    for i = 1:size(mixtures,1)
        if any(mixtures(i,1) == [1,6,7,10,11,19,20]) 
            bacnames{i,1} = bacNamegroups{1};
            if any(mixtures(i,3) == [1,6,7,10,11,19,20])
                bacnames{i,2} = bacNamegroups{1};
                bacgroups(i) = 1;
            elseif any(mixtures(i,3) ==[8,9,12:18,22])
                bacnames{i,2} = bacNamegroups{2};
                bacgroups(i) = 2;
            elseif any(mixtures(i,3) == [2,3,4,5,21])
                bacnames{i,2} = bacNamegroups{3};
                bacgroups(i) = 5;
            end 
        elseif any(mixtures(i,1) == [8,9,12:18,22]) 
             bacnames{i,1} = bacNamegroups{2};
            if any(mixtures(i,3) == [8,9,12:18,22])
                bacnames{i,2} = bacNamegroups{2};
                bacgroups(i)=4;
            elseif any(mixtures(i,1) == [2,3,4,5,21])
                bacnames{i,2} = bacNamegroups{3};
                bacgroups(i)=8;
            end 
        elseif any(mixtures(i,1) == [2,3,4,5,21]) && any(mixtures(i,1) == [2,3,4,5,21])
            bacnames{i,1} = bacNamegroups{3};
            bacnames{i,2} = bacNamegroups{3};
            bacgroups(i) = 9;
        else 
            bacnames{i,1} ='0';
        end 
    end 
end 