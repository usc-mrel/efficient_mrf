%% MRF pattern matching script.
clear; clc; close all;

addpath(genpath('../nih'));

data_dir = '../../data/volunteer/pulseq_mrf_20240329';
folder_name = {'nufft_3d', 'subspace_3d', 'maxgirf_3d'};

Ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
Nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];
opt_dict = '';
opt_data = '';
for opt_indx = [1 3 5 8 11 14] % Phantom, 1:14
    dict_folder = sprintf('../../data/dictionaries/disc/pulseq_readout_experiments/dict_ni%d_nt%d%s', Ni(opt_indx), Nt(opt_indx), opt_dict);
    
    tstart = tic;
    fprintf('Loading dictionaries with B1+ correction ...\n');
    ii = 1;
    for rf = 0.5:0.05:1.5
        dict_name = sprintf('fisp_mrf_3d_ni%d_nt%d%s_rf%g.mat', Ni(opt_indx), Nt(opt_indx), opt_dict, rf);
        dict_name = fullfile(dict_folder, dict_name);
        load(dict_name);
        if rf == 0.5
            [nt, nte, nr] = size(dict);
            dict_all = complex(zeros(nt, nr, length(0.5:0.05:1.5)));
        end
        dict_all(:, :, ii) = squeeze(dict(:, 2, :));
        ii = ii + 1;
    end
    fprintf('Done! Time cost: %.2f sec.\n', toc(tstart));
    
    for slice_loc = 20:40
        mat_names{1} = fullfile(data_dir, 'nufft_3d', sprintf('pulseq_mrf_slice%d__fa75_cycle2_ni%d_nt%d%s.mat', slice_loc, Ni(opt_indx), Nt(opt_indx), opt_data));
        mat_names{2} = fullfile(data_dir, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice%d_B1__fa75_cycle2_ni%d_nt%d%s.mat', slice_loc, Ni(opt_indx), Nt(opt_indx), opt_data));
        mat_names{3} = fullfile(data_dir, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice%d_B0__fa75_cycle2_ni%d_nt%d%s.mat', slice_loc, Ni(opt_indx), Nt(opt_indx), opt_data));
        
        N = [256 256];
        
        %% Load B1 map.
        b1_range = 0.5:0.05:1.5;
        if exist(fullfile(data_dir, 'b1mapping/b1map.mat'), 'file')
            fprintf('Loading B1 map ...\n');
            load(fullfile(data_dir, 'b1mapping/b1map.mat'));
        else
            fprintf('Rounding B1 map to the nearest 0.05 and convert B1 map to an index map ...\n');
            [b1map, b1_indx] = convert_dicom_to_b1(fullfile(data_dir, 'b1mapping'), b1_range.', 0);
        end
        b1_indx = rot90(flip(b1_indx, 1), -1);
        fprintf('Done! Time cost: %.2f sec.\n', toc(tstart));
        
        for ii = 1:length(mat_names)
            mat_name = mat_names{ii};
            load(mat_name);
            
            fprintf('Simple pattern matching for %s ...\n', mat_name);
            
            nblocks = 1;
            t1map = zeros([N, 1]);
            t2map = t1map;
            m0map = t2map;
            matches = m0map;
            r(:,4:5) = 0;
            
            switch folder_name{ii}
                case 'nufft_3d'
                    recon_im = recon_im;
                case 'subspace_3d'
                    recon_im = mrf_image_subspace;
                case 'maxgirf_3d'
                    recon_im = mrf_image_maxgirf;
            end
            
            for jj = 1:size(dict_all, 3) % Do pattern matching with all dictionaries
                t_dict = tic; fprintf('Performing pattern matching using dictionary with B1+=%0.2f.\n', b1_range(jj));
                indx = (b1_indx(:, :, slice_loc) == jj);
                [t1, t2, ~, m0, match] = patternmatch(squeeze(recon_im), indx, r, 0, dict_all(:, :, jj), nblocks);
                t1map = t1map + t1;
                t2map = t2map + t2;
                m0map = m0map + mat2gray(abs(m0));
                matches = matches + match;
                fprintf('Done! Time cost: %.2f sec.\n', toc(t_dict));
            end
            
            % Save the results.
            fprintf('Saving the results ...\n');
            save(mat_name, 't1map', 't2map', 'm0map', 'matches', '-append');
            fprintf('Done! (%.2f sec)\n', toc(tstart));
            
            clear recon_im mrf_image_subspace mrf_image_maxgirf;
        end
    end
end