clear; clc; close all

%%% Tight subplot path.
addpath(genpath('~/Documents/MATLAB/mrel/util/tight_subplot'));
addpath('~/Documents/MATLAB/mrel/3d_mrf/pulseq_mrf_recon');

ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];

%%% Load subspace and MaxGIRF reconstructed images and T1maps.
data_path = '/Users/zhibozhu/Documents/mrel_data/mrf/volunteer/readouts/pulseq_mrf_20240401';
subspace_path = dir(fullfile(data_path, './subspace_3d'));
maxgirf_path = dir(fullfile(data_path, './maxgirf_3d'));

tr = [7 8.1 9.1 10.2 11.1 12.4 13.6 15.1 16.1 17.3 18.8 20.6 23 26.1];
ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];

load('/Users/zhibozhu/Documents/mrel_data/mrf/T1cm.mat');
load('/Users/zhibozhu/Documents/mrel_data/mrf/T2cm.mat');

%% Reconstructed images.
load(fullfile(data_path, 'mask_slice25.mat'));
reorient = @(x) flip(rot90(x, -1), 2);
text_color = [255 204 0] ./ 255;

fS5 = figure('Position', [200 160 925 605], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

xrange = 116:145;
yrange = 151:200;
window_x = [116:145 145*ones(1, 50) 145:-1:116 116*ones(1, 50)];
window_y = [151*ones(1, 30) 151:200, 200*ones(1, 30) 200:-1:151];
inset_x = [ones(1, 50)-0.5 0.5:30.5 30.5*ones(1, 50) 30.5:-1:0.5];
inset_y = [0.5:50.5 50.5*ones(1, 30) 50.5:-1:0.5 ones(1, 30)-0.5];

ncol = 0;
for ii = [14 8 3]
    load(fullfile(data_path, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice25_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_subspace');
    load(fullfile(data_path, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice25_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        'mrf_image_maxgirf');

    nr_time_frame = round(100 * 7 / tr(15 - ii));
    mrf_image_subspace = abs(mrf_image_subspace(:, :, nr_time_frame));
    mrf_image_maxgirf = abs(mrf_image_maxgirf(:, :, nr_time_frame));

    ncol = ncol + 1;
    axes(ha(1, ncol));
    imagesc(reorient(mask .* mrf_image_subspace), [0 4]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [0 255 128] / 255);
    text(60, 70, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 13, 'FontWeight', 'bold');
    if ncol == 1
        text(50, 128, 'Gridding+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    if ii == 14 && ncol == 1
        ax1 = gca;
    else
        linkaxes([ax1, gca]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0.5  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange) .* mrf_image_subspace(xrange, yrange)), [0 4]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [0 255 128] / 255);

    axes(ha(2, ncol));
    imagesc(reorient(mask .* mrf_image_maxgirf), [0 4]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [0 255 128] / 255);
    linkaxes([ax1, gca]);
    if ncol == 1
        text(50, 128, 'MaxGIRF+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end
    if ncol == 3
        xlim([41 220]);
        ylim([51 230]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange) .* mrf_image_maxgirf(xrange, yrange)), [0 4]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [0 255 128] / 255);

end
colormap(gray);

%% T1 maps.
fS5_t1 = figure('Position', [200 160 1060 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

ncol = 0;
for ii = [14 8 3]
    load(fullfile(data_path, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice25_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't1map_v2');

    ncol = ncol + 1;
    axes(ha(1, ncol));
    imagesc(reorient(mask .* t1map_v2), [0 1.5e3]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [0 255 128] / 255);
    text(60, 70, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 13, 'FontWeight', 'bold');
    if ncol == 1
        text(50, 128, 'Gridding+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    if ii == 14 && ncol == 1
        ax1 = gca;
    else
        linkaxes([ax1, gca]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0.5  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange) .* t1map_v2(xrange, yrange)), [0 1.5e3]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [0 255 128] / 255);

    if ncol == 3
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 12;
        c.Label.String = '[ms]';
        c.Label.Position = [0.5 1650 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.9550 0.6700 0.0100 0.3000];
    end

    load(fullfile(data_path, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice25_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't1map_v2');
    axes(ha(2, ncol));
    imagesc(reorient(mask .* t1map_v2), [0 1.5e3]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [0 255 128] / 255);
    linkaxes([ax1, gca]);
    if ncol == 1
        text(50, 128, 'MaxGIRF+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end
    if ncol == 3
        xlim([41 220]);
        ylim([51 230]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange) .* t1map_v2(xrange, yrange)), [0 1.5e3]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [0 255 128] / 255);
end
colormap(T1colormap);

%% T2 maps.
fS5_t2 = figure('Position', [200 160 1060 625], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

ncol = 0;
for ii = [14 8 3]
    load(fullfile(data_path, 'subspace_3d', sprintf('pulseq_mrf_subspace_slice25_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't2map_v2');

    ncol = ncol + 1;
    axes(ha(1, ncol));
    imagesc(reorient(mask .* t2map_v2), [0 2e2]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [0 255 128] / 255);
    text(60, 70, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 13, 'FontWeight', 'bold');
    if ncol == 1
        text(50, 128, 'Gridding+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    if ii == 14 && ncol == 1
        ax1 = gca;
    else
        linkaxes([ax1, gca]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0.5  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange) .* t2map_v2(xrange, yrange)), [0 2e2]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [0 255 128] / 255);

    if ncol == 3
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 12;
        c.Label.String = '[ms]';
        c.Label.Position = [0.5 220 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.9550 0.6700 0.0100 0.3000];
    end

    load(fullfile(data_path, 'maxgirf_3d', sprintf('pulseq_mrf_maxgirf_slice25_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), ...
        't2map_v2');
    axes(ha(2, ncol));
    imagesc(reorient(mask .* t2map_v2), [0 2e2]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [0 255 128] / 255);
    linkaxes([ax1, gca]);
    if ncol == 1
        text(50, 128, 'MaxGIRF+Subspace', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end
    if ncol == 3
        xlim([41 220]);
        ylim([51 230]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange) .* t2map_v2(xrange, yrange)), [0 2e2]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [0 255 128] / 255);
end
colormap(T2colormap);

%% Save figures.
saveas(fS5, 'figureS5.fig');
saveas(fS5, 'figureS5.png');
saveas(fS5_t1, 'figureS5_t1.fig');
saveas(fS5_t1, 'figureS5_t1.png');
saveas(fS5_t2, 'figureS5_t2.fig');
saveas(fS5_t2, 'figureS5_t2.png');


