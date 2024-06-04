%% This is also Figure S7.
clear; clc; close all

%%% Tight subplot path.
addpath(genpath('../src/utils/tight_subplot'));
addpath('../src/pulseq_mrf_recon');

ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];

%%% Load subspace and MaxGIRF reconstructed images and T1maps.
data_path = '../data/volunteer/pulseq_mrf_20240329';
subspace_path = dir(fullfile(data_path, './subspace_3d'));
maxgirf_path = dir(fullfile(data_path, './maxgirf_3d'));

tr = [7 8.1 9.1 10.2 11.1 12.4 13.6 15.1 16.1 17.3 18.8 20.6 23 26.1];
ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];

load('../T1cm.mat');
load('../T2cm.mat');

%% Slice 36.
load(fullfile(data_path, 'mask_slice36.mat'));
reorient = @(x) flip(rot90(x, -1), 2);
text_color = [255 204 0] ./ 255;

f8_p1 = figure('Position', [200 160 925 605], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

xrange = 116:145;
yrange = 151:200;
window_x = [116:145 145*ones(1, 50) 145:-1:116 116*ones(1, 50)];
window_y = [151*ones(1, 30) 151:200, 200*ones(1, 30) 200:-1:151];
inset_x = [ones(1, 50)-0.5 0.5:30.5 30.5*ones(1, 50) 30.5:-1:0.5];
inset_y = [0.5:50.5 50.5*ones(1, 30) 50.5:-1:0.5 ones(1, 30)-0.5];

load(fullfile(data_path, 'mrf_vol1_analysis_slice36.mat'));
ncol = 0;
for ii = [14 8 3]
    switch ii
        case 14
            t1map_sd = t1map_maxgirf_std(:, :, 1);
            t2map_sd = t2map_maxgirf_std(:, :, 1);
        case 8
            t1map_sd = t1map_maxgirf_std(:, :, 3);
            t2map_sd = t2map_maxgirf_std(:, :, 3);
        case 3
            t1map_sd = t1map_maxgirf_std(:, :, 5);
            t2map_sd = t2map_maxgirf_std(:, :, 5);
    end
    ncol = ncol + 1;
    axes(ha(1, ncol));
    imagesc(reorient(mask) .* t1map_sd, [0 50]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [255 0 0] / 255);
    text(60, 70, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 18, 'FontWeight', 'bold');
    if ncol == 1
        text(50, 128, 'T_1 standard deviation', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    if ii == 14 && ncol == 1
        ax1 = gca;
    else
        linkaxes([ax1, gca]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0.5  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange)) .* t1map_sd(yrange, xrange), [0 50]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [255 0 0] / 255);

    if ncol == 3
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 16;
        c.Label.String = '[ms]';
        c.Label.Position = [0 55 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.9550 0.6700 0.0100 0.3000];
    end

    axes(ha(2, ncol));
    imagesc(reorient(mask) .* t2map_sd, [0 10]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [255 0 0] / 255);
    linkaxes([ax1, gca]);
    if ncol == 1
        text(50, 128, 'T_2 standard deviation', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end
    if ncol == 3
        xlim([41 220]);
        ylim([51 230]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange)) .* t2map_sd(yrange, xrange), [0 10]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [255 0 0] / 255);

    if ncol == 3
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 16;
        c.Label.String = '[ms]';
        c.Label.Position = [0 11 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.9550 0.1700 0.0100 0.3000];
    end
end
cm = colormap(jet);
cm(1, :) = [0 0 0];
colormap(cm);

%% Slice 25.
load(fullfile(data_path, 'mask_slice25.mat'));
reorient = @(x) flip(rot90(x, -1), 2);
text_color = [255 204 0] ./ 255;

fS7_v2 = figure('Position', [200 160 925 605], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 3, [0 0], [0 0], [0.03 0.03]);
ha = (reshape(ha, [3 2])).';

xrange = 116:145;
yrange = 151:200;
window_x = [116:145 145*ones(1, 50) 145:-1:116 116*ones(1, 50)];
window_y = [151*ones(1, 30) 151:200, 200*ones(1, 30) 200:-1:151];
inset_x = [ones(1, 50)-0.5 0.5:30.5 30.5*ones(1, 50) 30.5:-1:0.5];
inset_y = [0.5:50.5 50.5*ones(1, 30) 50.5:-1:0.5 ones(1, 30)-0.5];

load(fullfile(data_path, 'mrf_vol1_analysis_slice25.mat'));
ncol = 0;
for ii = [14 8 3]
    switch ii
        case 14
            t1map_sd = t1map_maxgirf_std(:, :, 1);
            t2map_sd = t2map_maxgirf_std(:, :, 1);
        case 8
            t1map_sd = t1map_maxgirf_std(:, :, 3);
            t2map_sd = t2map_maxgirf_std(:, :, 3);
        case 3
            t1map_sd = t1map_maxgirf_std(:, :, 5);
            t2map_sd = t2map_maxgirf_std(:, :, 5);
    end
    ncol = ncol + 1;
    axes(ha(1, ncol));
    imagesc(reorient(mask) .* t1map_sd, [0 50]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [255 0 0] / 255);
    text(60, 70, sprintf('%g ms', ro(15 - ii)), 'Color', text_color, 'FontSize', 16, 'FontWeight', 'bold');
    if ncol == 1
        text(50, 128, 'T_1 standard deviation', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end

    if ii == 14 && ncol == 1
        ax1 = gca;
    else
        linkaxes([ax1, gca]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0.5  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange)) .* t1map_sd(yrange, xrange), [0 50]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [255 0 0] / 255);

    if ncol == 3
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 16;
        c.Label.String = '[ms]';
        c.Label.Position = [0.5 55 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.9550 0.6700 0.0100 0.3000];
    end

    axes(ha(2, ncol));
    imagesc(reorient(mask) .* t2map_sd, [0 10]);
    axis off image;
    hold on;
    plot(window_x, window_y, '--', 'LineWidth', 2, 'Color', [255 0 0] / 255);
    linkaxes([ax1, gca]);
    if ncol == 1
        text(50, 128, 'T_2 standard deviation', 'FontSize', 18, 'FontWeight', 'bold', 'Color', text_color, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    end
    if ncol == 3
        xlim([41 220]);
        ylim([51 230]);
    end

    axes('Position', [(0.18 + (ncol - 1) * 0.3133) 0  0.2 0.2]);
    imagesc(reorient(mask(xrange, yrange)) .* t2map_sd(yrange, xrange), [0 10]);
    axis off image;
    hold on;
    plot(inset_x, inset_y, '-', 'LineWidth', 2, 'Color', [255 0 0] / 255);

    if ncol == 3
        pos = get(gca, 'Position');
        c = colorbar;
        set(gca, 'Position', pos);
        c.FontSize = 16;
        c.Label.String = '[ms]';
        c.Label.Position = [0.5 11 0];
        c.Label.Rotation = 0;
        c.Color = text_color;
        c.Position = [0.9550 0.1700 0.0100 0.3000];
    end
end
cm = colormap(jet);
cm(1, :) = [0 0 0];
colormap(cm);

%% Save figures.
saveas(f8_p1, 'figure8.fig');
saveas(f8_p1, 'figure8.png');
saveas(fS7_v2, 'figureS7_v2.fig');
saveas(fS7_v2, 'figureS7_v2.png');


