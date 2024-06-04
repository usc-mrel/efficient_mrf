clear; clc; close all

%%% Tight subplot path.
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

fS3 = figure('Position', [200 160 1060 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(4, 7, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [7 4])).';

for ii = 14:-1:1
    load(fullfile(data_path, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice36_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_subspace');
    load(fullfile(data_path, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_maxgirf');

    nrow = fix((14 - ii) / 7) + 1;
    ncol = rem((14 - ii), 7) + 1;
    if nrow == 2
        nrow = 3;
    end

    nr_time_frame = round(100 * 7 / tr(15 - ii));
    mrf_image_subspace = abs(mrf_image_subspace(:, :, nr_time_frame));
    mrf_image_maxgirf = abs(mrf_image_maxgirf(:, :, nr_time_frame));

    axes(ha(nrow, ncol));
    imagesc(reorient(mask .* mrf_image_subspace), [0 8]);
    axis off image;
%     hold on;
%     plot(window_x, window_y, 'r--', 'LineWidth', 2);
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 13, 'FontWeight', 'bold');
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'Gridding+Subspace', 'FontSize', 14, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

%     axes('Position', [(0.08 + (ncol - 1) * 0.1343) (0.75 - (nrow - 1) * 0.25)  0.1 0.1]);
%     imagesc(reorient(mask(xrange, yrange) .* mrf_image_subspace(xrange, yrange)), [0 8]);
%     axis off image;
%     hold on;
%     plot(inset_x, inset_y, 'r-', 'LineWidth', 2);

    axes(ha(nrow + 1, ncol));
    imagesc(reorient(mask .* mrf_image_maxgirf), [0 8]);
    axis off image;
%     hold on;
%     plot(window_x, window_y, 'r--', 'LineWidth', 2);
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'MaxGIRF+Subspace', 'FontSize', 14, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

%     axes('Position', [(0.08 + (ncol - 1) * 0.1343) (0.75 - nrow * 0.25)  0.1 0.1]);
%     imagesc(reorient(mask(xrange, yrange) .* mrf_image_maxgirf(xrange, yrange)), [0 8]);
%     axis off image;
%     hold on;
%     plot(inset_x, inset_y, 'r-', 'LineWidth', 2);
end
colormap(gray);

%% Edge maps.
fS4 = figure('Position', [200 160 1235 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(4, 7, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [7 4])).';

xrange = 146:165;
yrange = 82:101;
window_x = [146:165 165*ones(1, 20) 165:-1:146 146*ones(1, 20)];
window_y = [82*ones(1, 20) 82:101, 101*ones(1, 20) 101:-1:82];
inset_x = [ones(1, 20)-0.5 0.5:20.5 20.5*ones(1, 20) 20.5:-1:0.5];
inset_y = [0.5:20.5 20.5*ones(1, 20) 20.5:-1:0.5 ones(1, 20)-0.5];

for ii = 14:-1:1
    nrow = fix((14 - ii) / 7) + 1;
    ncol = rem((14 - ii), 7) + 1;
    if nrow == 2
        nrow = 3;
    end

    load(fullfile(data_path, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice36_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't1map');
    axes(ha(nrow, ncol));
    imagesc(reorient(mask .* t1map), [0 3e3]);
    axis off image;
    hold on;
    plot(window_x, window_y, 'r--', 'LineWidth', 2);
    text(50, 55, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 13, 'FontWeight', 'bold');
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'Gridding+Subspace', 'FontSize', 14, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    axes('Position', [(0.08 + (ncol - 1) * 0.1343) (0.75 - (nrow - 1) * 0.25)  0.1 0.1]);
    imagesc(reorient(mask(xrange, yrange) .* t1map(xrange, yrange)), [0 3e3]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, 'r-', 'LineWidth', 2);

    if nrow == 1 && ncol == 7
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 12;
        c.Label.String = '[ms]';
        c.Label.Position = [0.5 3300 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.965 0.75 0.01 0.23];
    end

    load(fullfile(data_path, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't1map');
    axes(ha(nrow + 1, ncol));
    imagesc(reorient(mask .* t1map), [0 3e3]);
    axis off image;
    hold on;
    plot(window_x, window_y, 'r--', 'LineWidth', 2);
    zoom(gca, 1.5);
    if ncol == 1
        text(30, 128, 'MaxGIRF+Subspace', 'FontSize', 14, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    axes('Position', [(0.08 + (ncol - 1) * 0.1343) (0.75 - nrow * 0.25)  0.1 0.1]);
    imagesc(reorient(mask(xrange, yrange) .* t1map(xrange, yrange)), [0 3e3]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, 'r-', 'LineWidth', 2);
end
colormap(T1colormap);

%% Save figures.
saveas(fS3, 'figureS3.fig');
saveas(fS3, 'figureS3.png');
saveas(fS4, 'figureS4.fig');
saveas(fS4, 'figureS4.png');

