clear; clc; close all

%%% Tight subplot path.
addpath(genpath('~/Documents/MATLAB/mrel/util/tight_subplot'));
% addpath('~/Documents/MATLAB/mrel/3d_mrf/pulseq_mrf_recon');

ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];

%%% Load spatial coordinates.
load('concomitant_field_ni30.mat', 'p');

tr = [7 8.1 9.1 10.2 11.1 12.4 13.6 15.1 16.1 17.3 18.8 20.6 23 26.1];
ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];

%% Figure 1, spiral arms and concomitant field induced phases.
gamma = 4257.59 * (1e4 * 2 * pi); % [rad/sec/T]
Nl = 19;
B0 = [0.55 1.5];
dz = [0 0.05 0.1 0.15];

clear phase;
f1 = figure('Position', [260 500 910 590/2], 'Color', 'w', 'InvertHardcopy', 'off');
text_color = [0.85 0.85 0.85];
[ha, pos] = tight_subplot(1, 3, [0.1 0.1], 2*[0.1 0.08], [0.08 0.05]);
ha = (reshape(ha, [3 1])).';

%%% Ni = 48.
ii = 0;
for Ni_indx = [14 8 3]
    ii = ii + 1;

    % 0.55T, Aera
    load(sprintf('traj_ni%d_nt%d.mat', ni(Ni_indx), nt(Ni_indx)));
    krmax = 2 * pi / (0.3 / 256) /2;
    k_rcs = kx_GIRF(:, 1) * 2 * krmax;
    k_rcs(:, :, 2) = -ky_GIRF(:, 1) * 2 * krmax;
    k_rcs(:, :, 3) = 0;
    k_rcs = permute(k_rcs, [3 1 2]);
    R_rcs2gcs = [0 1 0; 1 0 0; 0 0 1];
    k_gcs = R_rcs2gcs * k_rcs;
    g_gcs = diff(k_gcs, [], 2) / (4257.59e4 * 2 * pi * 1e-3 * 2.5e-6); % [mT/m]
    R_gcs2dcs = [0 1 0; -1 0 0; 0 0 -1];
    g_dcs = R_gcs2dcs * g_gcs * 1e-1;
    k1 = calculate_concomitant_field_coefficients(permute(squeeze(g_dcs(1, :)), [2 1]), permute(squeeze(g_dcs(2, :)), [2 1]), permute(squeeze(g_dcs(3, :)), [2 1]), Nl, B0(1), gamma, dt);
    k3 = calculate_concomitant_field_coefficients(permute(squeeze(g_dcs(1, :)), [2 1]), permute(squeeze(g_dcs(2, :)), [2 1]), permute(squeeze(g_dcs(3, :)), [2 1]), Nl, B0(2), gamma, dt);

    % 0.55T, Freemax
    load(sprintf('traj_ni%d_nt%d_freemax.mat', ni(Ni_indx), nt(Ni_indx)));
    krmax = 2 * pi / (0.3 / 256) /2;
    k_rcs = kx_GIRF(:, 1) * 2 * krmax;
    k_rcs(:, :, 2) = -ky_GIRF(:, 1) * 2 * krmax;
    k_rcs(:, :, 3) = 0;
    k_rcs = permute(k_rcs, [3 1 2]);
    R_rcs2gcs = [0 1 0; 1 0 0; 0 0 1];
    k_gcs = R_rcs2gcs * k_rcs;
    g_gcs = diff(k_gcs, [], 2) / (4257.59e4 * 2 * pi * 1e-3 * 2.5e-6); % [mT/m]
    R_gcs2dcs = [0 1 0; -1 0 0; 0 0 -1];
    g_dcs = R_gcs2dcs * g_gcs * 1e-1;
    k2 = calculate_concomitant_field_coefficients(permute(squeeze(g_dcs(1, :)), [2 1]), permute(squeeze(g_dcs(2, :)), [2 1]), permute(squeeze(g_dcs(3, :)), [2 1]), Nl, B0(1), gamma, dt);

    p(:, 3) = dz(3);
    p(:, 6) = p(:, 3).^2;
    p(:, 12) = p(:, 3).^3;

    temp = k1(:, 4:end) * p(:, 4:end).';
    phase1 = mean(temp, 2) / (2 * pi);
    temp = k2(:, 4:end) * p(:, 4:end).';
    phase2 = mean(temp, 2) / (2 * pi);
    temp = k3(:, 4:end) * p(:, 4:end).';
    phase3 = mean(temp, 2) / (2 * pi);

    axes(ha(1, ii));
    l1 = plot(dt * (0:size(phase1, 1)) * 1e3, [0; phase1], 'LineWidth', 2, 'Color', [0.92 0.694 0.125]);
    hold on;
    l2 = plot(dt * (0:size(phase2, 1)) * 1e3, [0; phase2], 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
    l3 = plot(dt * (0:size(phase3, 1)) * 1e3, [0; phase3], 'LineWidth', 2, 'Color', [0 0.447 0.741]);
    
    xlabel('Time [ms]', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase [cycles]', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([0 32]);
    ylim([0 8]);
    legend('0.55 T, Aera', '0.55 T, Free.MAX', '1.5 T, Aera', 'Location', 'northwest');
    title([sprintf('#Interleaves = %g, ', ni(Ni_indx)) '|\Deltaz| = 10 cm'], 'FontSize', 18);

    clear phase1 phase2 phase3;
end

%% Save figures.
% saveas(f1, 'figure1.fig');
saveas(f1, 'figure1_extra.png');


