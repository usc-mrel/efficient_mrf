clear; clc; close all;

load('../T1cm.mat');
load('../T2cm.mat');

ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];
reorient = @(x) x.';
text_color = [255 204 0] / 255;
f = figure('Position', [200 160 800 800], 'Color', 'k', 'InvertHardcopy', 'off');
[ha, pos] = tight_subplot(2, 2, [0 0], [0 0], [0.03 0.05]);
ha = (reshape(ha, [2 2])).';

load('../data/volunteer/pulseq_mrf_20240329/mrf_vol1_3d.mat');
load('../data/volunteer/pulseq_mrf_20240329/mask_3d.mat');

t1map_subspace = t1map_subspace .* mask;
t2map_subspace = t2map_subspace .* mask;
t1map_maxgirf = t1map_maxgirf .* mask;
t2map_maxgirf = t2map_maxgirf .* mask;

yrange = 36:225;
xrange = 51:210;
for ii = 1:33
    axes(ha(1, 1));
    ax1 = gca;
    imagesc(reorient(t1map_subspace(xrange, yrange, ii)), [0 1.5e3]);
    axis off image;
    text(-10, 10, ['\Delta z = ' sprintf('%g mm', (ii+7-24.5)*5)], 'Color', text_color, 'FontSize', 18, 'FontWeight', 'bold')
    text(-10, 20, 'RO = 16.5 ms', 'Color', text_color, 'FontSize', 18, 'FontWeight', 'bold');

    axes(ha(1, 2));
    ax2 = gca;
    imagesc(reorient(t1map_maxgirf(xrange, yrange, ii)), [0 1.5e3]);
    axis off image;

    pos = get(gca, 'Position');
    c1 = colorbar;
    set(gca, 'Position', pos);
    c1.FontSize = 12;
    c1.Label.String = '[ms]';
    c1.Label.Position = [0.5 1650 0];
    c1.Label.Rotation = 0;
    c1.Color = text_color;
    c1.Position = [0.95 0.6700 0.0100 0.3000];

    axes(ha(2, 1));
    ax3 = gca;
    imagesc(reorient(t2map_subspace(xrange, yrange, ii)), [0 2e2]);
    axis off image;
    
    axes(ha(2, 2));
    ax4 = gca;
    imagesc(reorient(t2map_maxgirf(xrange, yrange, ii)), [0 2e2]);
    axis off image;

    pos = get(gca, 'Position');
    c2 = colorbar;
    set(gca, 'Position', pos);
    c2.FontSize = 12;
    c2.Label.String = '[ms]';
    c2.Label.Position = [0.5 220 0];
    c2.Label.Rotation = 0;
    c2.Color = text_color;
    c2.Position = [0.950 0.1700 0.0100 0.3000];

    linkaxes([ax1 ax2 ax3 ax4]);
    % zoom(ax1, 1.35);

    colormap(ax1, T1colormap);
    colormap(ax2, T1colormap);
    colormap(ax3, T2colormap);
    colormap(ax4, T2colormap);

    im = frame2im(getframe(f));
    [A, map] = rgb2ind(im,256);
    if ii == 1
        imwrite(A, map, 'movieM3.gif', 'gif', 'LoopCount', Inf, 'DelayTime', 1);
    else
        imwrite(A, map, 'movieM3.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 1);
    end
end