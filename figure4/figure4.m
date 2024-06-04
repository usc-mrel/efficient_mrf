clear; clc; close all

%%% Tight subplot path.
addpath(genpath('../src/utils/tight_subplot'));

data_path = '../data/phantom/pulseq_mrf_20231113';

%% Ni = 48;
f4abc = figure('Position', [260 500 940 280], 'Color', 'k', 'InvertHardcopy', 'off');
text_color = [255 204 0] ./ 255;
[ha, pos] = tight_subplot(1, 3, [0 0.05], [0 0.05], [0.05 0.05]);
ha = (reshape(ha, [3 1])).';

load(fullfile(data_path, 'subspace_3d_iso/pulseq_mrf_subspace_slice24_B1__fa75_cycle2_ni48_nt1036.mat'));
mrf_image_subspace_iso = abs(mrf_image_subspace(:, :, 518));
% mrf_image_subspace_iso = mrf_image_subspace_iso / max(mrf_image_subspace_iso(:));
axes(ha(1, 1));
imagesc(mrf_image_subspace_iso.', [0 16]);
axis off square;
colormap gray;
title('(A) No corrections, \Deltaz=0mm, RO=2.9ms', 'FontSize', 14, 'Color', text_color);
ax1 = gca;

load(fullfile(data_path, 'subspace_3d/pulseq_mrf_subspace_slice24_B1__fa75_cycle2_ni48_nt1036.mat'));
mrf_image_subspace_h75mm = abs(mrf_image_subspace(:, :, 518));
% mrf_image_subspace_h75mm = mrf_image_subspace_h75mm / max(mrf_image_subspace_h75mm(:));
axes(ha(1, 2));
imagesc(mrf_image_subspace_h75mm.', [0 16]);
axis off square;
colormap gray;
title('(B) No corrections, \Deltaz=75mm, RO=2.9ms', 'FontSize', 14, 'Color', text_color);
ax2 = gca;

load(fullfile(data_path, 'maxgirf_3d/pulseq_mrf_maxgirf_slice24_B0__fa75_cycle2_ni48_nt1036.mat'));
mrf_image_maxgirf_h75mm = abs(mrf_image_maxgirf(:, :, 518));
% mrf_image_maxgirf_h75mm = mrf_image_maxgirf_h75mm / max(mrf_image_maxgirf_h75mm(:));
axes(ha(1, 3));
imagesc(mrf_image_maxgirf_h75mm.', [0 16]);
axis off square;
colormap gray;
title('(C) MaxGIRF correction, \Deltaz=75mm, RO=2.9ms', 'FontSize', 14, 'Color', text_color);
ax3 = gca;

linkaxes([ax1 ax2 ax3]);
zoom(ax1, 2);

%% Ni = 14;
f4def = figure('Position', [260 500 940 280], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(1, 3, [0 0.05], [0 0.05], [0.05 0.05]);
ha = (reshape(ha, [3 1])).';

load(fullfile(data_path, 'subspace_3d_iso/pulseq_mrf_subspace_slice24_B1__fa75_cycle2_ni14_nt535.mat'));
mrf_image_subspace_iso = abs(mrf_image_subspace(:, :, 267));
% mrf_image_subspace_iso = mrf_image_subspace_iso / max(mrf_image_subspace_iso(:));
axes(ha(1, 1));
imagesc(mrf_image_subspace_iso.', [0 16]);
axis off square;
colormap gray;
title('(D) No corrections, \Deltaz=0mm, RO=9.5ms', 'FontSize', 14, 'Color', text_color);
ax1 = gca;

load(fullfile(data_path, 'subspace_3d/pulseq_mrf_subspace_slice24_B1__fa75_cycle2_ni14_nt535.mat'));
mrf_image_subspace_h75mm = abs(mrf_image_subspace(:, :, 267));
% mrf_image_subspace_h75mm = mrf_image_subspace_h75mm / max(mrf_image_subspace_h75mm(:));
axes(ha(1, 2));
imagesc(mrf_image_subspace_h75mm.', [0 16]);
axis off square;
colormap gray;
title('(E) No corrections, \Deltaz=75mm, RO=9.5ms', 'FontSize', 14, 'Color', text_color);
ax2 = gca;

load(fullfile(data_path, 'maxgirf_3d/pulseq_mrf_maxgirf_slice24_B0__fa75_cycle2_ni14_nt535.mat'));
mrf_image_maxgirf_h75mm = abs(mrf_image_maxgirf(:, :, 267));
% mrf_image_maxgirf_h75mm = mrf_image_maxgirf_h75mm / max(mrf_image_maxgirf_h75mm(:));
axes(ha(1, 3));
imagesc(mrf_image_maxgirf_h75mm.', [0 16]);
axis off square;
colormap gray;
title('(F) MaxGIRF correction, \Deltaz=75mm, RO=9.5ms', 'FontSize', 14, 'Color', text_color);
ax3 = gca;

linkaxes([ax1 ax2 ax3]);
zoom(ax1, 2);

%% Save figures;
saveas(f4abc, 'figure4abc.fig');
saveas(f4abc, 'figure4abc.png');
saveas(f4def, 'figure4def.fig');
saveas(f4def, 'figure4def.png');


