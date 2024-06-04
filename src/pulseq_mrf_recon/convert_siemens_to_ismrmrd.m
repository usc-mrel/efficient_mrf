clear; clc; close all;

data_dir = '../../data/volunteer/pulseq_mrf_20240401';
cur_dir = pwd;

cd(data_dir);

twix_list = dir('meas*.dat');

for ii = 1:length(twix_list)
    twix_name = twix_list(ii).name;
    twix_name = twix_name(1:end-4);
    
    linux_command_noise = sprintf('siemens_to_ismrmrd -f %s.dat -z 1 -o noise_%s.h5', twix_name, twix_name);
    linux_command_data = sprintf('siemens_to_ismrmrd -f %s.dat -z 2 -o %s.h5', twix_name, twix_name);

    system(linux_command_noise);
    system(linux_command_data);
end

cd(cur_dir);