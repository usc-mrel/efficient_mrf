clear; clc; close all

%%% Tight subplot path.
addpath(genpath('../src/utils/tight_subplot'));
addpath('../src/pulseq_mrf_recon');
addpath(genpath('../src/utils/tight_subplot'));
addpath('../src/pulseq_mrf_recon');

ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];

%%% Load subspace and MaxGIRF reconstructed images and T1maps.
data_path = '../data/phantom/pulseq_mrf_20240306';
subspace_path = dir(fullfile(data_path, './subspace_3d'));
maxgirf_path = dir(fullfile(data_path, './maxgirf_3d'));

tr = [7 8.1 9.1 10.2 11.1 12.4 13.6 15.1 16.1 17.3 18.8 20.6 23 26.1];
ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];

load('../T1cm.mat');

%% Reconstructed images.
load(fullfile(data_path, 'mask.mat'));
reorient = @(x) flip(rot90(x, -1), 2);
text_color = [255 204 0] ./ 255;

f5 = figure('Position', [200 160 1060 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

xrange = 146:165;
yrange = 82:101;
window_x = [146:165 165*ones(1, 20) 165:-1:146 146*ones(1, 20)];
window_y = [82*ones(1, 20) 82:101, 101*ones(1, 20) 101:-1:82];
inset_x = [ones(1, 20)-0.5 0.5:20.5 20.5*ones(1, 20) 20.5:-1:0.5];
inset_y = [0.5:20.5 20.5*ones(1, 20) 20.5:-1:0.5 ones(1, 20)-0.5];

ncol = 0;
for ii = [14 8 3]
    load(fullfile(data_path, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice36_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_subspace');
    load(fullfile(data_path, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_maxgirf');

    nr_time_frame = round(100 * 7 / tr(15 - ii));
    mrf_image_subspace = abs(mrf_image_subspace(:, :, nr_time_frame));
    mrf_image_maxgirf = abs(mrf_image_maxgirf(:, :, nr_time_frame));

    ncol = ncol + 1;
    axes(ha(1, ncol));
    imagesc(reorient(mask .* mrf_image_subspace), [0 8]);
    axis off image;
    hold on;
    plot(window_x, window_y, 'r--', 'LineWidth', 2);
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 16, 'FontWeight', 'bold');
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'Gridding+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    axes('Position', [(0.22 + (ncol - 1) * 0.3133) 0.5  0.15 0.15]);
    imagesc(reorient(mask(xrange, yrange) .* mrf_image_subspace(xrange, yrange)), [0 8]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, 'r-', 'LineWidth', 2);

    axes(ha(2, ncol));
    imagesc(reorient(mask .* mrf_image_maxgirf), [0 8]);
    axis off image;
    hold on;
    plot(window_x, window_y, 'r--', 'LineWidth', 2);
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'MaxGIRF+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    axes('Position', [(0.22 + (ncol - 1) * 0.3133) 0  0.15 0.15]);
    imagesc(reorient(mask(xrange, yrange) .* mrf_image_maxgirf(xrange, yrange)), [0 8]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, 'r-', 'LineWidth', 2);
end
colormap(gray);

%% T1 maps.
f6 = figure('Position', [200 160 1060 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

ncol = 0;
for ii = [14 8 3]
    load(fullfile(data_path, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice36_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't1map');

    ncol = ncol + 1;
    axes(ha(1, ncol));
    imagesc(reorient(mask .* t1map), [0 3e3]);
    axis off image;
    hold on;
    plot(window_x, window_y, 'r--', 'LineWidth', 2);
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 16, 'FontWeight', 'bold');
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'Gridding+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    axes('Position', [(0.22 + (ncol - 1) * 0.3133) 0.5  0.15 0.15]);
    imagesc(reorient(mask(xrange, yrange) .* t1map(xrange, yrange)), [0 3e3]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, 'r-', 'LineWidth', 2);

    if ncol == 3
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 12;
        c.Label.String = '[ms]';
        c.Label.Position = [0.5 3300 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.9550 0.6700 0.0100 0.3000];
    end

    load(fullfile(data_path, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't1map');
    axes(ha(2, ncol));
    imagesc(reorient(mask .* t1map), [0 3e3]);
    axis off image;
    hold on;
    plot(window_x, window_y, 'r--', 'LineWidth', 2);
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'MaxGIRF+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    axes('Position', [(0.22 + (ncol - 1) * 0.3133) 0  0.15 0.15]);
    imagesc(reorient(mask(xrange, yrange) .* t1map(xrange, yrange)), [0 3e3]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, 'r-', 'LineWidth', 2);
end
colormap(T1colormap);

%% Reconstructed images, movie.
movieM2 = figure('Position', [200 160 1235 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(1, 2, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [2 1])).';

for ii = 14:-1:1
    load(fullfile(subspace_path(ii).folder, sprintf('pulseq_mrf_subspace_slice36_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_subspace');
    load(fullfile(maxgirf_path(ii).folder, sprintf('pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_maxgirf');

    nr_time_frame = round(100 * 7 / tr(15 - ii));
    mrf_image_subspace = abs(mrf_image_subspace(:, :, nr_time_frame));
    mrf_image_maxgirf = abs(mrf_image_maxgirf(:, :, nr_time_frame));

    axes(ha(1, 1));
    ax1 = gca;
    imagesc(reorient(mask .* mrf_image_subspace), [0 8]);
    axis off image;
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 24, 'FontWeight', 'bold');

    axes(ha(1, 2));
    ax2 = gca;
    imagesc(reorient(mask .* mrf_image_maxgirf), [0 8]);
    axis off image;

    linkaxes([ax1 ax2]);
    zoom(ax1, 1.5);
    colormap('gray');

    im = frame2im(getframe(f3));
    [A, map] = rgb2ind(im,256);
    if ii == 14
        imwrite(A, map, 'movieM2.gif', 'gif', 'LoopCount', Inf, 'DelayTime', 0.5);
    else
        imwrite(A, map, 'movieM2.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end

%% Save figures.
saveas(f5, 'figure5.fig');
saveas(f5, 'figure5.png');
saveas(f6, 'figure6.fig');
saveas(f6, 'figure6.png');


