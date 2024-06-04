function dict_pulseq(fn_MRF_seq_params, fn_MRF_seq_params_mat, output_dir)
% clear; clc; close all;
%
% fn_MRF_seq_params = 'mrf_params_ni48_nt1036.csv';
% fn_MRF_seq_params_mat = 'fisp_mrf_3d_seq_params_ni48_nt1036.mat';
load(fn_MRF_seq_params_mat, 'minTR', 'minTE', 'nr_time_frames', 'nr_interleaves');

% dictionary params
% T1 = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000]; % [ms]
% T2 = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500]; %[ms]
T1 = [10:10:2000 2050:50:3000 3100:100:3500];
T2 = [2:2:100 105:5:300 310:10:500 550:50:800 900:100:2000];
df = 50;

% setup dictionary parameters
data_struct = importdata(fn_MRF_seq_params);
B1 = 0.5:0.05:1.5;
TI = 20.6 / 1e3; % [s]
FA = data_struct.data(:, 1) * 75 * pi / 180; % [rad]
phi = data_struct.data(:, 2) * pi / 180; % [rad]
TR = minTR + data_struct.data(:, 3); % [s]
TE = [minTE minTE];
% TE = minTE + data_struct.data(:, 4); % [s]
waiting_delay = 2; % [s]

%
output_dir = sprintf('%s/dict_%s', output_dir, fn_MRF_seq_params_mat(82:end-4));
if ~exist(output_dir)
    mkdir(output_dir);
end
dephasing = 4;
inv = true;
r = paramTable(T1, T2, df);

cnt = size(r, 1);
p = gcp('nocreate');
if isempty(p)
    p = gcp;
end
poolsize = p.NumWorkers;
blocksize = ceil(cnt/poolsize);

for i = 1:length(B1)
    b1 = B1(i);
    dictr = zeros(nr_time_frames * length(TE), blocksize, poolsize);
    dicti = dictr;
    parfor j = 1:poolsize
        r_index = (j - 1)* blocksize + 1:j * blocksize;
        r_index(r_index > cnt) = [];
        dr = zeros(nr_time_frames * length(TE), blocksize);
        di = dr;
        [dr(:, 1:length(r_index)), di(:, 1:length(r_index)), ~] = dictionary_FISP_withRelaxDecay(...
            r(r_index,:), b1 * FA, TR, TE, phi, TI, dephasing * pi, 250, nr_time_frames, inv, waiting_delay);
        dictr(:,:,j) = dr;
        dicti(:,:,j) = di;
        
    end
    dict = complex(dictr, dicti);
    dict = reshape(dict, [nr_time_frames, length(TE), blocksize*poolsize]);
    dict = dict(:,:,1:cnt);
    
%     dict = dict(:,2:end,:); % throw away all the previous repeats, and the first echo signal
%     dict = circshift(dict,1,2); % shift the last echo to the first
    fn_dict_output = fullfile(output_dir, sprintf('fisp_mrf_3d_%s_df50_rf%g.mat', fn_MRF_seq_params(70:end-4), b1));
    
    % TO-DO: Align the dictionary with readout phase.
    dict = dict .* repmat(exp(-1i*phi(1:nr_time_frames)), [1 2 cnt]) .* 1i;
    
    save(fn_dict_output, 'dict', 'r');
end

% delete(p);