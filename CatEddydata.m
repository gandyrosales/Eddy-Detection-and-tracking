%% catEddydata before Sorting (CCSEddySort)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS
% in Rosales-Quintana et al -- 2024, september.
%
% To concatenate the Detection data before using CCSEddySort.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathin = '/Path_to_detected_outputs/';
info = dir(fullfile(pathin,'name_detected_eddies*.mat'));

datas = '';
for i = 1:1:size(info,1)
    disp(info(i).name);
    fn = fullfile(pathin,info(i).name);
    load(fn)
    datas = cat(2,datas,data);
end

% Here we Sort the concatenated data 
Eddies = CCSEddySort(datas);
% save('Eddies.mat', '-v7.3','-nocompression');

