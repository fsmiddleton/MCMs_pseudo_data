%% Batch processing 
% FS Middleton 10/11/2022
%Reference colors 
%Diana (2022). Color blind friendly colormap (https://www.mathworks.com/matlabcentral/fileexchange/46802-color-blind-friendly-colormap), MATLAB Central File Exchange. Retrieved November 18, 2022.
%% Adding files and starting 
clc
clear

c = parcluster('HPC1');
Temps = [298.15 288.15 313.15 318.15 323.15];
for i = 1:length(Temps)
    file{i} = strcat('heUNIFACforT=',num2str(Temps(i)),'.mat'); 
end 
%% 3way     
numworkers = 31;
job = batch(c, @Completion3wayParallelScript, 10, {}, 'Pool', numworkers, 'AttachedFiles', file); 
%% 2way 
c = parcluster('HPC1');
job = batch(c, @Completion2wayParallelScript, 11, {}, 'Pool', 31, 'AttachedFiles', file);
%% 2way unfolding
c = parcluster('HPC1');
job = batch(c, @Completion2wayParallelScriptUnfold, 11, {}, 'Pool', 31, 'AttachedFiles', file);
%% Postprocess 2way or 3way 

job_output=job168_output; % change the output file here 

filenamesave = job_output{1,1};
Xm_boot = job_output{1,2};
Xm_boot2 = job_output{1,3};
fns = job_output{1,5};
mse_LOOCV = job_output{1,4};
wmse_LOOCV = job_output{1,6};
X = job_output{1,7};
Xs = job_output{1,8};
conc_interval = job_output{1,9};
filename = job_output{1,10};
%Factors =  job_output{1,11};
save(filenamesave)
load(filenamesave)

%% Retrieve jobs 
c = parcluster('HPC1');
[pending queued running completed] = findJob(c);

i=1; %number of the completed job 
jobfetched = findJob(c,'ID',completed(i,1).ID);
job_output = fetchOutputs(jobfetched);


%%
c = parcluster('local');
c = parcluster('HPC1');
delete(c.Jobs);




