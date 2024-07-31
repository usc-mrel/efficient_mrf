% The GPU version is fast but has memory risk related to the MaxGIRF part.
% A100 preferred.
% Zhibo.

%% Clean slate
clear; clc; close all;

%% Set source directories
src_directory        = './src';
utils_directory      = '../utils';
export_fig_directory = '../utils/export_fig';
pulseq_directory     = '../../pulseq/pulseq';
ismrmrd_directory    = '../src/utils/ismrmrd/matlab';

%% Add source directories to search path
addpath(genpath(src_directory));
addpath(utils_directory);
addpath(genpath(export_fig_directory));
addpath(genpath(pulseq_directory));
addpath(genpath(ismrmrd_directory));

%% Set file paths
Nr_interleaves      = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
Nt                  = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];
opt_data            = '';
opt_dict            = '';
opt_traj            = '';
for opt_indx            = [1 3 5 8 11 14] % Phantom, 1:14
    seq_path            = sprintf('../../pulseq/lowfield_mrf/different_ro_fa_fm/fisp_mrf_3d_fa75_cycle2_ni%d_nt%d.seq', Nr_interleaves(opt_indx), Nt(opt_indx));
    traj_dir            = sprintf('../../pulseq/pulseq/lowfield_mrf/different_ro_fa_fm/traj_vdspiral_ni%d_nt%d%s.h5', Nr_interleaves(opt_indx), Nt(opt_indx), opt_traj);
    data_dir            = '../../data/volunteer/pulseq_mrf_20240329';
    dat_name            = sprintf('meas_MID*_FID*_fisp_mrf_3d_fa75_cycle2_ni%d_nt%d%s.dat', Nr_interleaves(opt_indx), Nt(opt_indx), opt_data);
    dict_path           = sprintf('../../data/dictionaries/disc/pulseq_20230705/dict_ni%d_nt%d%s', Nr_interleaves(opt_indx), Nt(opt_indx), opt_dict);

    %% Read a .seq file
    start_time = tic;
    tstart = tic; fprintf('Reading a .seq file: %s... ', seq_path);
    seq = mr.Sequence;
    seq.read(seq_path);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Get sequence parameters
    base_resolution = seq.getDefinition('BaseResolution');
    discard_post    = seq.getDefinition('DiscardPost');
    discard_pre     = seq.getDefinition('DiscardPre');
    real_dwell_time = seq.getDefinition('RealDwellTime'); % [sec]

    %% Define traj directory
    traj_dset = ismrmrd.Dataset(traj_dir, 'dataset');
    traj_header = ismrmrd.xml.deserialize(traj_dset.readxml);

    %% Export FoV, matrix size and saved trajectories
    encoded_fov(1) = traj_header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
    encoded_fov(2) = traj_header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
    encoded_fov(3) = traj_header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

    recon_fov(1) = traj_header.encoding.reconSpace.fieldOfView_mm.x * 1e-3; % [m] RO
    recon_fov(2) = traj_header.encoding.reconSpace.fieldOfView_mm.y * 1e-3; % [m] PE
    recon_fov(3) = traj_header.encoding.reconSpace.fieldOfView_mm.z * 1e-3; % [m] SL

    Nkx = traj_header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space, Cartesian
    Nky = traj_header.encoding.encodedSpace.matrixSize.y; % number of phase encodes in k-space, Cartesian
    Nkz = traj_header.encoding.encodedSpace.matrixSize.z; % number of slice encodes in k-space, Cartesian
    Nx  = traj_header.encoding.reconSpace.matrixSize.x;   % number of samples in image-space (RO)
    Ny  = traj_header.encoding.reconSpace.matrixSize.y;   % number of samples in image-space (PE)
    Nz  = traj_header.encoding.reconSpace.matrixSize.z;   % number of samples in image-space (SL)

    encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [m]
    recon_resolution = recon_fov ./ [Nx Ny Nz]; % [m]

    raw_traj = traj_dset.readAcquisition();
    for ii = 1:raw_traj.head.scan_counter(end)
        if ii == 1
            k_gcs_GRT = zeros([size(raw_traj.traj{ii}), raw_traj.head.scan_counter(end)]);
        end
        k_gcs_GRT(:, :, ii) = raw_traj.traj{ii};
    end

    %% Define data directory
    % twix_list = dir(fullfile(data_dir, 'dat_h75mm', dat_name)); % ACR MRI phantom
    twix_list = dir(fullfile(data_dir, dat_name)); % Others

    %% Batch recon
    verbose = 0;
    for ii = 1:length(twix_list)
        twix_name           = twix_list(ii).name;
        twix_name           = twix_name(1:end-4); % Trim the suffix, Zhibo.
        siemens_twix_path   = sprintf('%s/%s.dat', data_dir, twix_name);
        ismrmrd_data_path   = sprintf('%s/%s.h5', data_dir, twix_name);
        ismrmrd_noise_path  = sprintf('%s/noise_%s.h5', data_dir, twix_name);

        % For data mimicking Free.Max, downsample the readout by 2 to avoid memory issues.
        if contains(twix_name, 'freemax')
            readout_os_factor = 2;
        else
            readout_os_factor = 1;
        end

        % Read k-space data (ISMRMRD format)
        start_time = tic;
        tic;
        fprintf('Reading an ISMRMRD file: %s... ', ismrmrd_data_path);
        if exist(ismrmrd_data_path, 'file')
            dset = ismrmrd.Dataset(ismrmrd_data_path, 'dataset');
            fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
        else
            error('File %s does not exist. Please generate it.' , ismrmrd_data_path);
        end

        % Get imaging parameters from the XML header
        header = ismrmrd.xml.deserialize(dset.readxml);

        Nc  = header.acquisitionSystemInformation.receiverChannels;

        %--------------------------------------------------------------------------
        % measurement information
        %--------------------------------------------------------------------------
        patient_position = header.measurementInformation.patientPosition;

        % Parse the ISMRMRD header
        tic; fprintf('Parsing the ISMRMRD header... ');
        raw_data = dset.readAcquisition(); % read all the acquisitions
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

        grad_samples       = double(max(raw_traj.head.number_of_samples));
        adc_samples        = double(max(raw_data.head.number_of_samples));
        center_sample      = double(max(raw_data.head.center_sample));
        nr_channels        = double(max(raw_data.head.active_channels));
        nr_phase_encodings = double(max(raw_data.head.idx.kspace_encode_step_1)) + 1; % nr_interleaves for spiral imaging
        nr_slice_encodings = double(max(raw_data.head.idx.kspace_encode_step_2)) + 1;
        % nr_averages        = double(max(raw_data.head.idx.average)) + 1;
        % nr_slices          = double(max(raw_data.head.idx.slice)) + 1;
        nr_slices          = 48; % Hardcoded.
        nr_contrasts       = double(max(raw_data.head.idx.contrast)) + 1;
        nr_phases          = double(max(raw_data.head.idx.phase)) + 1;
        % nr_repetitions     = double(max(raw_data.head.idx.repetition)) + 1;
        nr_sets            = double(max(raw_data.head.idx.set)) + 1;
        nr_segments        = double(max(raw_data.head.idx.segment)) + 1;
        nr_samples         = adc_samples - discard_pre - discard_post;
        nr_averages        = seq.getDefinition('Averages');
        nr_repetitions     = seq.getDefinition('Repetitions');

        % Get the dimensionality of the trajectory vector (0 means no trajectory)
        trajectory_dimensions = double(max(raw_data.head.trajectory_dimensions));

        % Get the dwell time in [sec]
        % dt = double(max(raw_data.head.sample_time_us)) * 1e-6; % [usec] * [sec/1e-6 usec] => [sec]
        dt  = 2.5e-6; % Hardcoded.

        % Calculate the readout duration [sec]
        T   = adc_samples * dt; % readout duration [sec]

        % Update ISMRMRD encoding parameters
        % Calculate reconstruction parameters
        Nk  = Nkx; % number of readout samples
        N1  = Nk;
        N2  = Nky;
        N3  = Nkz;
        N   = N1 * N2 * N3;

        % Display ISMRMRD header
        if verbose
            fprintf('========================= ISMRMRD header ========================\n');
            fprintf('encoded_fov        = %8.4f %8.4f %8.4f\n', encoded_fov(1), encoded_fov(2), encoded_fov(3));
            fprintf('Nkx Nky Nkz        = %d      %d        %d\n', Nkx, Nky, Nkz);
            fprintf('encoded_resolution = %8.4f %8.4f %8.4f\n', encoded_resolution(1), encoded_resolution(2), encoded_resolution(3));
            fprintf('-----------------------------------------------------------------\n');
            fprintf('recon_fov          = %8.4f %8.4f %8.4f\n', recon_fov(1), recon_fov(2), recon_fov(3));
            fprintf('Nx Ny Nz           = %d      %d        %d\n', Nx, Ny, Nz);
            fprintf('recon_resolution   = %8.4f %8.4f %8.4f\n', recon_resolution(1), recon_resolution(2), recon_resolution(3));
            fprintf('-----------------------------------------------------------------\n');
            fprintf('trajectory         = %s\n', header.encoding.trajectory);
            fprintf('number_of_samples  = %d\n', adc_samples);
            fprintf('discard_pre        = %d\n', discard_pre);
            fprintf('discard_post       = %d\n', discard_post);
            fprintf('center_sample      = %d\n', center_sample);
            fprintf('nr_channels        = %d\n', nr_channels);
            fprintf('nr_phase_encodings = %d\n', nr_phase_encodings);
            fprintf('nr_slice_encodings = %d\n', nr_slice_encodings);
            fprintf('nr_averages        = %d\n', nr_averages);
            fprintf('nr_slices          = %d\n', nr_slices);
            fprintf('nr_contrasts       = %d\n', nr_contrasts);
            fprintf('nr_phases          = %d\n', nr_phases);
            fprintf('nr_repetitions     = %d\n', nr_repetitions);
            fprintf('nr_sets            = %d\n', nr_sets);
            fprintf('nr_segments        = %d\n', nr_segments);
            fprintf('dt                 = %5.2f [usec]\n', dt * 1e6);
            fprintf('readout duration   = %5.2f [msec]\n', T * 1e3);
            fprintf('=================================================================\n');
        end

        % Calculate the receiver noise matrix
        [Psi,inv_L] = calculate_receiver_noise_matrix(ismrmrd_noise_path);

        % Read a Siemens .dat file
        if exist(siemens_twix_path, 'file')
            fprintf('Reading a Siemens .dat file: %s\n', siemens_twix_path);
            twix = mapVBVD(siemens_twix_path);
            if length(twix) > 1
                twix = twix{end};
            end
        end

        % Get a slice normal vector from Siemens TWIX format
        % dNormalSag: Sagittal component of a slice normal vector (in the PCS)
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dSag')
            dNormalSag  = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
        else
            dNormalSag  = 0;
        end

        % dNormalCor: Coronal component of a slice normal vector (in the PCS)
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dCor')
            dNormalCor  = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
        else
            dNormalCor  = 0;
        end

        % dNormalTra: Transverse component of a slice normal vector (in the PCS)
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dTra')
            dNormalTra  = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
        else
            dNormalTra  = 0;
        end

        % dRotAngle: Slice rotation angle ("swap Fre/Pha")
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'dInPlaneRot')
            dRotAngle   = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot; % [rad]
        else
            dRotAngle   = 0; % [rad]
        end

        % Determine the main orientation of an imaging stack
        main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);

        % Calculate a scaling matrix [m]
        scaling_matrix  = diag(encoded_resolution); % [m]

        % Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
        R_rcs2gcs       = [0    1    0 ; % [PE]   [0 1 0] * [r]
            1    0    0 ; % [RO] = [1 0 0] * [c]
            0    0    1]; % [SL]   [0 0 1] * [s]

        % Calculate a rotation matrix from the GCS to the PCS
        [R_gcs2pcs, phase_sign, read_sign] = ...
            siemens_calculate_matrix_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

        % Calculate a rotation matrix from the PCS to the DCS
        R_pcs2dcs       = siemens_calculate_matrix_pcs_to_dcs(patient_position);

        % Calculate a rotation matrix from the GCS to the DCS
        R_gcs2dcs       = R_pcs2dcs * R_gcs2pcs;

        % Calculate a rotation matrix from the RCS to the DCS
        R_rcs2dcs       = R_pcs2dcs * R_gcs2pcs * R_rcs2gcs;

        % Define parameters for non-Cartesian reconstruction
        % Calculate the number of samples per echo
        Nk              = nr_samples;

        %--------------------------------------------------------------------------
        % Calculate the number of spiral interleaves
        %--------------------------------------------------------------------------
        Ni              = seq.getDefinition('Interleaves');

        % Calculate trajetories
        grad_raster_time = double(max(raw_data.head.sample_time_us)) * 1e-6;
        g_gcs_GRT       = diff(cat(2, zeros(3, 1, raw_traj.head.scan_counter(end)), k_gcs_GRT), 1, 2) ...
            / (seq.sys.gamma * 2 * pi) / grad_raster_time * 1e3; % [mT/m]

        % Interpolate nominal gradient waveforms in the GCS [mT/m] [PE,RO,SL]
        t_GRT           = ((0:grad_samples-1) + 0.5).' * grad_raster_time;
        t_ADC           = ((0:Nk-1) + 0.5).' * real_dwell_time;

        g_gcs_nominal   = zeros(3, Nk, Ni, 'single');
        for idx1 = 1:Ni
            g_gcs_nominal(1,:,idx1) = interp1(t_GRT, g_gcs_GRT(1,:,idx1), t_ADC, 'linear', 'extrap');
            g_gcs_nominal(2,:,idx1) = interp1(t_GRT, g_gcs_GRT(2,:,idx1), t_ADC, 'linear', 'extrap');
            g_gcs_nominal(3,:,idx1) = interp1(t_GRT, g_gcs_GRT(3,:,idx1), t_ADC, 'linear', 'extrap');
        end

        % Calculate GIRF-predicted gradient waveforms in the GCS [mT/m] [PE,RO,SL]
        tRR             = 0; % custom clock-shift
        sR.R            = R_gcs2dcs;
        sR.T            = header.acquisitionSystemInformation.systemFieldStrength_T;
        [~, g_gcs]      = apply_GIRF(permute(g_gcs_nominal, [2 3 1]), real_dwell_time, sR, tRR); % k:[cycle/cm] and g:[G/cm]
        g_gcs           = permute(g_gcs, [3 1 2]); % Nk x Ni x 3 => 3 x Nk x Ni

        % Change the sign of GIRF-predicted gradient waveforms in the GCS [mT/m] [PE,RO,SL]
        g_gcs(1,:)      = phase_sign * g_gcs(1,:); % [mT/m] PE (gu)
        g_gcs(2,:)      = read_sign  * g_gcs(2,:); % [mT/m] RO (gv)

        % Additional sign flip if gradient waveforms are obtained from MRD format
        g_gcs(1,:,:)    = -g_gcs(1,:,:);

        % Calculate GIRF-predicted k-space trajectories in the GCS [rad/m] [PE,RO,SL]
        % Numerically integrate the coefficients
        % [Hz/T] * [2*pi rad/cycle] * [mT/m] * [T/1e3 mT] * [sec] => 2 * pi * 1e-3 [rad/m]
        k_gcs           = cumsum(seq.sys.gamma * g_gcs * real_dwell_time * (2 * pi * 1e-3), 2); % 3x Nk x Ni [rad/m]

        % Calculate GIRF-predicted k-space trajectories in the DCS [rad/m] [x,y,z]
        k_dcs = zeros(3, Nk, Ni, 'double');
        for i = 1:Ni
            k_dcs(:,:,i) = R_gcs2dcs * k_gcs(:,:,i); % 3 x Nk
        end

        % Calculate GIRF-predicted k-space trajectories in the RCS [rad/m] [R,C,S]
        k_rcs           = zeros(3, Nk, Ni, 'double');
        for i = 1:Ni
            k_rcs(:,:,i) = R_rcs2gcs.' * k_gcs(:,:,i); % Nk x 3
        end

        % Calculate the maximum k-space value [rad/m]
        krmax           = 2 * pi / (encoded_fov(1) / Nx) / 2; % [rad/m]

        % Trim k-space trajectories
        kxy             = squeeze(complex(k_rcs(1, :, :), k_rcs(2, :, :)));
        [kmax, idx]     = max(abs(kxy(:, 1)));
        k_rcs           = k_rcs(:, 1:idx, :);
        k_rcs           = k_rcs / (2 * krmax); % [-0.5, 0.5]

        % Calculate a density compensation function (1 x Nk x NiNs)
        tstart = tic;
        fprintf('Calculating a density compensation function ... ');
        g               = permute(squeeze(complex(g_gcs(1, 1:readout_os_factor:idx, 1), g_gcs(2, 1:readout_os_factor:idx, 1))), [2 1]);
        k               = permute(squeeze(complex(k_gcs(1, 1:readout_os_factor:idx, 1), k_gcs(2, 1:readout_os_factor:idx, 1))), [2 1]);
        w_full          = abs(repmat(g(:, 1), [Ni 1])) .* abs(sin(angle(repmat(g(:, 1), [Ni 1])) - angle(repmat(k(:, 1), [Ni 1]))));
        w_full          = w_full .* (1 / sum(w_full(:))); % NkNi x 1
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        % Prepare k-space data
        nr_time_frames  = nr_phase_encodings / (nr_repetitions * nr_averages * nr_slices);
        kspace          = complex(zeros(idx, nr_slices, nr_channels, nr_time_frames));
        kspace_discard  = complex(zeros(discard_pre, nr_slices, nr_channels, nr_time_frames));
        nr_interleaves  = mod(1:nr_time_frames, Ni);
        nr_interleaves(nr_interleaves == 0) = Ni;
        tstart          = tic;
        for nr_readouts = 1:nr_phase_encodings
            fprintf('Loading kspace data and performing noise pre-whitening (%5d/%5d) ...\n', nr_readouts, nr_phase_encodings);
            temp_data   = raw_data.data{1, nr_readouts};
            temp_data   = (inv_L * temp_data.').';

            discard_data = temp_data(1:discard_pre, :);
            discard_data = permute(discard_data, [1 3 2]);

            temp_data   = temp_data(discard_pre+1:end-discard_post, :);
            temp_data   = permute(temp_data, [1 3 2]);


            nr_time     = mod(nr_readouts, nr_time_frames);
            nr_block    = floor((nr_readouts-1) / nr_time_frames);

            % The following is for multi-repetition data only.
            %         if nr_repetitions > 1
            %             nr_interleaves = mod(nr_readouts - nr_block * nr_time_frames + nr_block, Ni);
            %         else
            %             nr_interleaves = mod(nr_readouts - nr_block * nr_time_frames, Ni);
            %         end

            % 3D data formatting.
            if nr_time ~= 0
                kspace(:, nr_block+1, :, nr_time) = temp_data(1:idx, :, :);
                kspace_discard(:, nr_block+1, :, nr_time) = discard_data;
            else
                kspace(:, nr_block+1, :, nr_time_frames) = temp_data(1:idx, :, :);
                kspace_discard(:, nr_block+1, :, nr_time_frames) = discard_data;
            end

            % 2D data formatting.
            %         if nr_interleaves ~= 0 && nr_time ~=0
            %             kspace(:, nr_interleaves, :, nr_time) = kspace(:, nr_interleaves, :, nr_time) + temp_data(1:idx, :, :);
            %         elseif nr_interleaves == 0 && nr_time ~= 0
            %             kspace(:, Ni, :, nr_time) = kspace(:, Ni, :, nr_time) + temp_data(1:idx, :, :);
            %         elseif nr_interleaves ~= 0 && nr_time == 0
            %             kspace(:, nr_interleaves, :, nr_time_frames) = kspace(:, nr_interleaves, :, nr_time_frames) + temp_data(1:idx, :, :);
            %         else
            %             kspace(:, Ni, :, nr_time_frames) = kspace(:, Ni, :, nr_time_frames) + temp_data(1:idx, :, :);
            %         end
        end
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time)');

        % Downsampling
        kspace          = kspace(1:readout_os_factor:end, :, :, :);
        k_rcs           = k_rcs(:, 1:readout_os_factor:end, :);
        Nk              = size(kspace, 1);

        % Averaging
        kspace          = kspace / nr_averages;

        % NUFFT recon
        run ../utils/mirt-master/setup;

        % Read GIRF corrected from the previous section.
        kx_GIRF         = -squeeze(k_rcs(1, :, :));
        ky_GIRF         = -squeeze(k_rcs(2, :, :));

        % Make NUFFT struct from MIRT toolbox.
        Nxy             = [Nx Ny];
        J               = [5 5];
        K               = 2 * Nxy;
        %     nshots_fullysampled = size(kx_GIRF, 2);
        %     tempindex       = 0:(nshots_fullysampled-1);
        %     trajindex       = zeros(length(tempindex), nshots_fullysampled);

        tempk           = [kx_GIRF(:) ky_GIRF(:)];
        nufft_st{1,1}   = nufft_init(2 * pi * tempk, Nxy, J, K, Nxy / 2, 'minmax:kb');
        DCF             = w_full;

        % GPU-NUFFT
        [~, device_idx] = gpuDeviceCount("available");
        gpuDevice(device_idx(3));

        nufft_st_device             = nufft_st;
        nufft_st_device{1, 1}.alpha = {gpuArray(nufft_st{1, 1}.alpha{1}), gpuArray(nufft_st{1, 1}.alpha{2})};
        nufft_st_device{1, 1}.beta  = {gpuArray(nufft_st{1, 1}.beta{1}), gpuArray(nufft_st{1, 1}.beta{2})};
        nufft_st_device{1, 1}.om    = gpuArray(nufft_st{1, 1}.om);
        nufft_st_device{1, 1}.sn    = gpuArray(nufft_st{1, 1}.sn);
        % nufft_st_device{1, 1}.p = gpuArray(nufft_st{1 ,1}.p); % This on cannot be transferred to GPU memory.
        DCF_device                  = gpuArray(DCF);


        % Transform 3d kspace into kx-ky-z.
        kspace          = ifftshift(ifft(kspace, [], 2), 2);

        %% NUFFT
        slice_loc       = [25 36];
        tstart          = tic;
        for nr_slice = slice_loc % 1:nr_slices if recon all slices.
            for nr_time = 1:nr_time_frames
                fprintf('Performing NUFFT reconstruction (%4d/%4d)(%2d/%2d)(%6.4f/%6.4f sec) ...\n', ...
                    nr_time, nr_time_frames, nr_slice, nr_slices, toc(tstart), toc(start_time));
                kspace_2d = complex(zeros(Nk, Ni, nr_channels));
                kspace_2d(:, nr_interleaves(nr_time), :) = kspace(:, nr_slice, :, nr_time);
                kspace_2d = reshape(kspace_2d, [], Nc);
                kspace_2d = gpuArray(kspace_2d);
                uncombined_img = nufft_adj_gpu(kspace_2d .* ...
                    repmat(DCF_device(:) * prod(nufft_st_device{1, 1}.Kd), [1 Nc]), ...
                    nufft_st_device{1, 1});

                if nr_slice == 1 && nr_time ==1
                    mrf_image_uncombined = complex(zeros(Nx, Ny, Nc, nr_time_frames));
                end

                mrf_image_uncombined(:, :, :, nr_time) = gather(uncombined_img);
            end

            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time)');

            clear uncombined_img

            t_avging    = mean(mrf_image_uncombined, 4);
            csm         = ismrm_estimate_csm_walsh(t_avging);
            clear tavging;
            csm_sq      = sum(csm .* conj(csm), 3);
            csm_sq(csm_sq < eps) = 1;

            recon_im    = complex(zeros(Nx, Ny, size(kspace, 4)));
            tsart       = tic;
            fprintf('Performing coil combination ...');
            for nr_time = 1:nr_time_frames
                recon_im(:, :, nr_time) = ((sqrt(numel(DCF(:)) / prod(K))) * ...
                    (sum(conj(csm) .* mrf_image_uncombined(:, :, :, nr_time), 3) ...
                    ./ csm_sq) ./ sqrt(prod(K)));
            end
            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

            % Create the output directories and save the recon images.
            if ~exist(fullfile(data_dir, 'csm'), 'dir')
                mkdir(fullfile(data_dir, 'csm'));
            end
            if ~exist(fullfile(data_dir, 'csm', sprintf('csm_ch%d_slice%d.mat', nr_channels, nr_slice)), 'file')
                save(fullfile(data_dir, 'csm', sprintf('csm_ch%d_slice%d.mat', nr_channels, nr_slice)), 'csm');
            end
            if ~exist(fullfile(data_dir, 'nufft_3d'), 'dir')
                mkdir(fullfile(data_dir, 'nufft_3d'));
            end
            save(fullfile(data_dir, 'nufft_3d', sprintf('pulseq_mrf_slice%d_%s.mat', nr_slice, twix_name(35:end))), 'recon_im', '-v7.3');
            clear recon_im mrf_image_uncombined;
        end

        %% Subspace learning.
        if ~exist(fullfile(dict_path, 'Phi.mat'), 'file')
            dict_names  = dir(fullfile(dict_path, '*rf*.mat'));
            for jj = 1:length(dict_names)
                load(fullfile(dict_names(jj).folder, dict_names(jj).name));

                % dict  = permute(output_dict.dict_norm, [2 1]); % This is ss-EPG approach for 2D MRF.
                dict    = permute(squeeze(dict(:, 2, :)), [2 1]); % This is for CWRU's approach.

                if jj == 1
                    dict_all = dict;
                else
                    dict_all = cat(1, dict_all, dict);
                end
            end

            Lt          = 25; % Specify the rank number.
            [U, S, V]   = svd(dict, 'econ');
            %         Phi         = S(1:Lt, 1:Lt) * (V(:, 1:Lt))';
            Phi         = V(:, 1:Lt)';
            save(fullfile(dict_path, 'Phi.mat'), 'Phi');
        else
            load(fullfile(dict_path, 'Phi.mat'));
            Lt = size(Phi, 1);
        end
        Phi_device      = gpuArray(Phi);

        % slice_loc = [15 20 24 25 29 34];
        for nr_slice = slice_loc
            % Load the CSM from the previous section.
            load(fullfile(data_dir, 'csm', sprintf('csm_ch%d_slice%d.mat', nr_channels, nr_slice)), 'csm');
            csm         = gpuArray(csm);

            % Form a zero-padded ksapce for the current slice.
            kspace_2d   = complex(zeros(size(kspace, 1), Ni, nr_channels, nr_time_frames));
            for nr_time = 1:nr_time_frames
                kspace_2d(:, nr_interleaves(nr_time), :, nr_time) = kspace(:, nr_slice, :, nr_time);
            end

            % Calculate b which is the RHS of the normal equation, AcH(Ac(U)) = AcH(d).
            Omega       = logical(abs(kspace_2d));
            kspace_us   = reshape(kspace_2d(Omega), size(kspace_2d, 1), 1, Nc, nr_time_frames);
            Omega       = gpuArray(reshape(Omega(:, :, 1, :), [], nr_time_frames));
            b           = calculate_b_subspace_gpu(Phi_device, kspace_us, Omega, csm, nufft_st_device, DCF);

            % Perform low-rank and subspace learning reconstruction.
            tol         = 1e-5;
            max_iterations = 100;
            start_time_lowrank_subspace = tic;
            fprintf('Performing low-rank and subspace model reconstruction...\n');
            E           = @(x,tr) encoding_non_cartesian_lowrank_subspace_gpu(x, csm, Omega, Phi_device, nufft_st_device, DCF_device, tr);
            [U_lowrank_subspace, flag, relres, iter, resvec] = lsqr(E, b, tol, max_iterations, [], [], []); % NLt x 1
            computation_time_lowrank_subspace = toc(start_time_lowrank_subspace);
            fprintf('Done! (%6.4f/%6.4f sec)\n', computation_time_lowrank_subspace, toc(start_time_lowrank_subspace));
            clear encoding_non_cartesian_lowrank_subspace_gpu;

            U_lowrank_subspace = gather(U_lowrank_subspace);
            mrf_image_subspace = reshape(U_lowrank_subspace, [], Lt) * Phi;
            mrf_image_subspace = reshape(mrf_image_subspace, [Nxy, nr_time_frames]);
            subspace_name = sprintf('pulseq_mrf_subspace_slice%d_B1_%s.mat', nr_slice, twix_name(35:end));
            if ~exist(fullfile(data_dir, 'subspace_3d'), 'dir')
                mkdir(fullfile(data_dir, 'subspace_3d'));
            end
            save(fullfile(data_dir, 'subspace_3d', subspace_name), 'U_lowrank_subspace', 'mrf_image_subspace', 'resvec');

            clear kspace_us b % Release memories
        end

        %% MaxGIRF
        % Get spatial coordinates in the DCS [x,y,z]
        r_dcs = zeros(N, 3, 'double');
        actual_slice_nr = 1;
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'sPosition')
            if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dSag')
                sag_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dSag; % [mm]
            else
                sag_offset_twix = 0; % [mm]
            end
            if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dCor')
                cor_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dCor; % [mm]
            else
                cor_offset_twix = 0; % [mm]
            end
            if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dTra')
                tra_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dTra; % [mm]
            else
                tra_offset_twix = 0; % [mm]
            end
        else
            sag_offset_twix = 0; % [mm]
            cor_offset_twix = 0; % [mm]
            tra_offset_twix = 0; % [mm]
        end
        pcs_offset      = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]

        % Get a slice offset in the PCS from ISMRMRD format
        sag_offset_ismrmrd = double(raw_data.head.position(1,actual_slice_nr)); % [mm]
        cor_offset_ismrmrd = double(raw_data.head.position(2,actual_slice_nr)); % [mm]
        tra_offset_ismrmrd = double(raw_data.head.position(3,actual_slice_nr)); % [mm]
        pcs_offset_ismrmrd = [sag_offset_ismrmrd; cor_offset_ismrmrd; tra_offset_ismrmrd] * 1e-3; % [mm] * [m/1e3mm] => [m]

        % Get a rotation matrix from the GCS to the PCS (ISMRMRD format)
        phase_dir       = double(raw_data.head.phase_dir(:,actual_slice_nr));
        read_dir        = double(raw_data.head.read_dir(:,actual_slice_nr));
        slice_dir       = double(raw_data.head.slice_dir(:,actual_slice_nr));
        R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

        % Calculate a slice offset in the DCS [m]
        dcs_offset      = R_pcs2dcs * pcs_offset; % 3 x 1

        % Calculate spatial coordinates in the RCS [m]
        [I1,I2,I3]      = ndgrid((1:N1).', (1:N2).', (1:N3).');
        r_rcs           = (scaling_matrix * cat(2, I1(:) - (floor(N1/2) + 1), I2(:) - (floor(N2/2) + 1), I3(:) - (floor(N3/2) + 1)).').'; % N x 3

        r_gcs           = (R_rcs2gcs * r_rcs.').'; % N x 3

        % Calculate spatial coordinates in the DCS [m]
        r_dcs(:,:)      = (repmat(dcs_offset, [1 N]) + R_pcs2dcs * R_gcs2pcs * r_gcs.').'; % N x 3

        g_dcs           = zeros(size(g_gcs, 2), 3, Ni, 'double');
        g_gcs           = permute(g_gcs, [2 1 3]) * 1e-1; % 1 [mT/m] -> 0.1 [G/cm]
        for i = 1:Ni
            g_dcs(:,:,i) = (R_gcs2dcs * g_gcs(:,:,i).').'; % Nk x 3
        end
        g_gcs           = permute(g_gcs, [2 1 3]) * 10; % 1 [G/cm] -> 10 [mT/m]
        g_dcs           = g_dcs(1:readout_os_factor:idx, :, :);

        gamma           = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]
        Nl              = 19; % Number of basis in the concomitant fields.
        B0              = 0.55; % Main field strength.

        % NUFFT struct per interleave.
        st              = cell(Ni,1);
        for i = 1:Ni
            tstart      = tic;
            fprintf('(%2d/%2d): Initializing structure for NUFFT per interleaf... ', i, Ni);
            om          = cat(2, kx_GIRF(:, i), ky_GIRF(:, i)); % Nk x 2
            %         Nd          = [Nx Ny];   % matrix size
            %         Jd          = [6 6];     % kernel size
            %         Kd          = Nd * 2;    % oversampled matrix size
            st{i}       = nufft_init(2 * pi * om, Nxy, J, K, Nxy/2, 'minmax:kb');

            st{i}.alpha = {gpuArray(st{i}.alpha{1}), gpuArray(st{i}.alpha{2})};
            st{i}.beta  = {gpuArray(st{i}.beta{1}), gpuArray(st{i}.beta{2})};
            st{i}.om    = gpuArray(st{i}.om);
            st{i}.sn    = gpuArray(st{i}.sn);

            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        end
        start_time      = tic;
        tstart          = tic;
        fprintf('Calculating the time courses of phase coefficients... ');
        k               = calculate_concomitant_field_coefficients(reshape(g_dcs(:,1,:), [size(g_dcs, 1) Ni]), reshape(g_dcs(:,2,:), [size(g_dcs, 1) Ni]), reshape(g_dcs(:,3,:), [size(g_dcs, 1) Ni]), Nl, B0, gamma, dt);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        % B0 map.
        fname_B0        = fullfile(data_dir, 'b0mapping/b0map.mat');
        % fullfile(data_dir, 'multi_gre_3d/nlinv_estimation_b0_3d/b0_gre_3d.mat');
        if ~isfile(fname_B0)
            load(fname_B0, 'b0map');
        else
            b0map       = zeros(Nx, Ny, nr_slices);
        end

        % Calculate the low-rank approximated higher order encoding operators.
        Lmax            = 12; % Rank number of the low rank approximation of the higher order encoding operators.
        os              = 5;
        t               = (0:readout_os_factor:idx-1)' * dt;
        b0map           = imresize3d(b0map, [Nx, Ny, nr_slices]);
        nr_recons       = nr_slices; % nr_slices
        static_B0_correction = 1;

        %     slice_loc = [15 20 24 25 29 34];
        for nr_slice = slice_loc
            r_dcs(:, 3) = abs(nr_slice - nr_slices / 2 - 0.5) * 0.005;
            %         r_dcs(:, 3) = 0; % Isocenter.
            %         r_dcs(:, 3) = 0.075; % H 75mm.
            tstart      = tic;
            printf('Calculating concomitant field basis functions... ');
            p           = calculate_concomitant_field_basis(r_dcs(:,1), r_dcs(:,2), r_dcs(:,3), Nl);
            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

            start_time_svd = tic;
            Ue          = complex(zeros(size(g_dcs, 1), Lmax, Ni, 'double', 'gpuArray'));
            Ve          = complex(zeros(Nx*Ny, Lmax, Ni, 'double', 'gpuArray'));
            s           = zeros(Lmax, Ni, 'double', 'gpuArray');
            for i = 1:Ni
                tstart = tic; fprintf('(%2d/%2d): Calculating randomized SVD (i=%2d/%2d)... ', nr_slice, nr_recons, i, Ni);
                [U_,S_,V_]  = calculate_rsvd_higher_order_encoding_matrix(k(:,4:end,i), p(:,4:end), Lmax, os, reshape(b0map(:, :, nr_slice), [Nx*Ny 1]), t, static_B0_correction); % TBD.
                Ue(:,:,i)   = U_(:,1:Lmax); % U: Nk x Lmax+os => Nk x Lmax
                Ve(:,:,i)   = V_(:,1:Lmax) * S_(1:Lmax,1:Lmax)'; % V: N x Lmax+os => N x Lmax
                s(:,i)      = diag(S_(1:Lmax,1:Lmax));
                fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
            end
            computation_time_svd = toc(start_time_svd);

            % Sanity checkL: This is equivalent ot subspace learning recon but much slower.
            % Ue              = ones(2680, 1, 12);
            % Ve              = ones(22500, 1, 12);

            % Load pre-calculated CSM.
            load(fullfile(data_dir, 'csm', sprintf('csm_ch%d_slice%d.mat', nr_channels, nr_slice)), 'csm');
            csm         = gpuArray(csm);

            kspace_2d   = complex(zeros(size(kspace, 1), Ni, nr_channels, nr_time_frames));
            for nr_time = 1:nr_time_frames
                kspace_2d(:, nr_interleaves(nr_time), :, nr_time) = kspace(:, nr_slice, :, nr_time);
            end
            Omega       = logical(abs(kspace_2d));
            Omega       = reshape(Omega(:, :, 1, :), [], nr_time_frames);
            Omega       = reshape(Omega, [], Ni, nr_time_frames);
            Omega       = permute(Omega, [1 3 2]); % Ns x Nt x Ni.
            Omega       = gpuArray(Omega);
            b           = calculate_b_maxgirf_subspace_gpu(Phi_device, Ue, Ve, kspace_2d, Omega, csm, st, reshape(DCF, [size(g_dcs, 1) Ni]));

            tol         = 1e-5;
            max_iterations = 50;
            start_time_maxgirf = tic;
            fprintf('(%d/%d): Performing CG-based MaxGIRF reconstruction...\n', slice_loc, nr_recons);
            E           = @(x,tr) encoding_lowrank_maxgirf_subspace_gpu(x, csm, Omega, Phi_device, Ue, Ve, reshape(DCF, [size(g_dcs, 1) Ni]), st, tr);
            [U_maxgirf, flag, relres, iter, resvec] = lsqr(E, b, tol, max_iterations, [], [], []);
            computation_time_maxgirf = toc(start_time_maxgirf);
            fprintf('done! (%6.4f/%6.4f sec)\n', computation_time_maxgirf, toc(start_time));
            clear encoding_lowrank_maxgirf_subspace_gpu;

            U_maxgirf = gather(U_maxgirf);
            mrf_image_maxgirf = reshape(U_maxgirf, [], Lt) * Phi;
            mrf_image_maxgirf = reshape(mrf_image_maxgirf, [Nxy nr_time_frames]);
            if ~static_B0_correction
                maxgirf_name = sprintf('pulseq_mrf_maxgirf_slice%d_%s_%s.mat', nr_slice, 'noB0', twix_name(35:end));
            else
                maxgirf_name = sprintf('pulseq_mrf_maxgirf_slice%d_%s_%s.mat', nr_slice, 'B0', twix_name(35:end));
            end
            if ~exist(fullfile(data_dir, 'maxgirf_3d'), 'dir')
                mkdir(fullfile(data_dir, 'maxgirf_3d'));
            end
            save(fullfile(data_dir, 'maxgirf_3d', maxgirf_name), 'U_maxgirf', 'mrf_image_maxgirf', 'resvec');
        end

        save(fullfile(data_dir, sprintf('traj_ni%d_nt%d%s.mat', Ni, Nt(opt_indx), opt_data)), 'kx_GIRF', 'ky_GIRF', 'dt');
        % Display the recon.
        verbose = 0;
        if verbose
            figure('Color', 'w');
            hold on;
            grid on; grid minor;
            set(gca, 'Box', 'On');
            X = reshape(r_dcs(:,1), [N1 N2 N3]);
            Y = reshape(r_dcs(:,2), [N1 N2 N3]);
            Z = reshape(r_dcs(:,3), [N1 N2 N3]);
            surf(X * 1e3, Y * 1e3, Z * 1e3, abs(mean(recon_im, 3)), 'EdgeColor', 'none');
            axis image;
            xlabel('x (R => L) [mm]');
            ylabel('y (P => A) [mm]');
            zlabel('z (H => F) [mm]');
            colormap(gray(256));
            set(gca, 'ZDir', 'reverse');
            if main_orientation == 2 % transversal
                view(0,90);
            elseif main_orientation == 1 % coronal
                view(0,0);
            elseif main_orientation == 0 % sagittal
                view(-90,0);
            end
            xlim([-150 150]);
            ylim([-150 150]);
            zlim([-150 150]);
            title(sprintf('PulSeq spiral MRF dephasing %s cycles', twix_name(48)));
            export_fig(fullfile(data_dir, sprintf(twix_name(24:end))), '-r300', '-tif');
        end
    end
end
