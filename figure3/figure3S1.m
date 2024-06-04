clear; clc; close all

%%% Tight subplot path.
addpath(genpath('../src/utils/tight_subplot'));
addpath('../src/pulseq_mrf_recon');

ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];

%%% Load subspace and MaxGIRF reconstructed images and T1maps.
data_path = '../data/phantom/pulseq_mrf_20231113';
subspace_path = dir(fullfile(data_path, './subspace_3d'));
maxgirf_path = dir(fullfile(data_path, './maxgirf_3d'));

tr = [7 8.1 9.1 10.2 11.1 12.4 13.6 15.1 16.1 17.3 18.8 20.6 23 26.1];
ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];

%% Reconstructed images.
load(fullfile(data_path, 'mask.mat'));
reorient = @(x) flip(rot90(x, -1), 2);
text_color = [255 204 0] ./ 255;

f3 = figure('Position', [200 160 1060 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

ncol = 0;
for ii = [14 8 3]
    load(fullfile(subspace_path(ii).folder, sprintf('pulseq_mrf_subspace_slice24_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_subspace');
    load(fullfile(maxgirf_path(ii).folder, sprintf('pulseq_mrf_maxgirf_slice24_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_maxgirf');

    nr_time_frame = round(100 * 7 / tr(15 - ii));
    mrf_image_subspace = abs(mrf_image_subspace(:, :, nr_time_frame));
    mrf_image_maxgirf = abs(mrf_image_maxgirf(:, :, nr_time_frame));

%     nrow = fix((14 - ii) / 7) + 1;
%     ncol = rem((14 - ii), 7) + 1;
%     if nrow == 2
%         nrow = 3;
%     end

    ncol = ncol + 1;
    axes(ha(1, ncol));
    ax1 = gca;
    imagesc(reorient(mask .* mrf_image_subspace), [0 16]);
    axis off image;
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 16, 'FontWeight', 'bold');
    zoom(ax1, 1.5);
    if ncol == 1
        text(30, 128, 'Gridding+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    if ii ~= 14
        linkaxes([ax1 ax2]);
    end

    axes(ha(2, ncol));
    ax2 = gca;
    imagesc(reorient(mask .* mrf_image_maxgirf), [0 16]);
    axis off image;
    zoom(ax2, 1.5);
    if ncol == 1
        text(30, 128, 'MaxGIRF+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    linkaxes([ax1 ax2]);
end
colormap('gray');

%% Edge maps.
f_not_used = figure('Position', [200 160 1060 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

ncol = 0;
for ii = [14 8 3]
    load(fullfile(subspace_path(ii).folder, sprintf('pulseq_mrf_subspace_slice24_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_subspace');
    load(fullfile(maxgirf_path(ii).folder, sprintf('pulseq_mrf_maxgirf_slice24_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_maxgirf');

    nr_time_frame = round(100 * 7 / tr(15 - ii));
    mrf_image_subspace = abs(mrf_image_subspace(:, :, nr_time_frame));
    mrf_image_maxgirf = abs(mrf_image_maxgirf(:, :, nr_time_frame));

%     nrow = fix((14 - ii) / 7) + 1;
%     ncol = rem((14 - ii), 7) + 1;
%     if nrow == 2
%         nrow = 3;
%     end

    ncol = ncol + 1;
    axes(ha(1, ncol));
    ax1 = gca;
    edge1 = mrf_image_subspace - imgaussfilt(mrf_image_subspace);
    imagesc(reorient(mask .* edge1), [-0.5 0.5]);
    axis off image;
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 16, 'FontWeight', 'bold');
    zoom(ax1, 1.5);
    if ncol == 1
        text(30, 128, 'NUFFT+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    if ii ~= 14
        linkaxes([ax1 ax2]);
    end

    axes(ha(2, ncol));
    ax2 = gca;
    edge2 = mrf_image_maxgirf - imgaussfilt(mrf_image_maxgirf);
    imagesc(reorient(mask .* edge2), [-0.5 0.5]);
    axis off image;
    zoom(ax2, 1.5);
    if ncol == 1
        text(30, 128, 'MaxGIRF+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    linkaxes([ax1 ax2]);
end
colormap('gray');

%% Reconstructed images, movie.
mM1 = figure('Position', [200 160 1235 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(1, 2, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [2 1])).';

for ii = 14:-1:1
    load(fullfile(subspace_path(ii).folder, sprintf('pulseq_mrf_subspace_slice24_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_subspace');
    load(fullfile(maxgirf_path(ii).folder, sprintf('pulseq_mrf_maxgirf_slice24_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_maxgirf');

    nr_time_frame = round(100 * 7 / tr(15 - ii));
    mrf_image_subspace = abs(mrf_image_subspace(:, :, nr_time_frame));
    mrf_image_maxgirf = abs(mrf_image_maxgirf(:, :, nr_time_frame));

    axes(ha(1, 1));
    ax1 = gca;
    imagesc(reorient(mask .* mrf_image_subspace), [0 16]);
    axis off image;
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 24, 'FontWeight', 'bold');

    axes(ha(1, 2));
    ax2 = gca;
    imagesc(reorient(mask .* mrf_image_maxgirf), [0 16]);
    axis off image;

    linkaxes([ax1 ax2]);
    zoom(ax1, 1.5);
    colormap('gray');

    im = frame2im(getframe(f3));
    [A, map] = rgb2ind(im,256);
    if ii == 14
        imwrite(A, map, 'movieM1.gif', 'gif', 'LoopCount', Inf, 'DelayTime', 0.5);
    else
        imwrite(A, map, 'movieM1.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end

%% Save figures.
saveas(f3, 'figure3.fig');
saveas(f3, 'figure3.png');
saveas(f_not_used, 'figure_not_used.fig');
saveas(f_not_used, 'figure_not_used.png');


