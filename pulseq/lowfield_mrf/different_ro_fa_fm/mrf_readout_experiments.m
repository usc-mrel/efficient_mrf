%% Clean slate
close all; clear; clc;

%% Set source directories
pulseq_directory = '../../pulseq';
thirdparty_directory = '../../thirdparty';
ismrmrd_directory    = '../../ismrmrd/matlab';

%% Add source directories to search path
addpath(genpath(pulseq_directory));
addpath(genpath(thirdparty_directory));
addpath(genpath(ismrmrd_directory));

%%
dirs = dir('./mrf_params_ni*_nt*.csv');
for ii = 1:length(dirs)
    fn_mrf_seq_path = fullfile(dirs(ii).folder, dirs(ii).name);
    temp = split(dirs(ii).name, '_');
    nr_interleaves = temp{3};
    nr_interleaves = str2double(nr_interleaves(3:end));
    make_mrf_sequence(nr_interleaves, fn_mrf_seq_path, 'freemax');
end

%%
dirs = dir('fisp_mrf_3d_seq_params_ni*_nt*.mat');
for ii = 1:length(dirs)
    load(fullfile(dirs(ii).folder, dirs(ii).name), 'minTR');
    disp(mean(minTR));
end
