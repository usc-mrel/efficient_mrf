clear; clc; close all;

load('../data/volunteer/pulseq_mrf_20240329/mrf_vol1_analysis_slice36.mat');
v = nrrdread('../data/volunteer/pulseq_mrf_20240329/seg_mprage_wm.nrrd');
v = double(v(:, :, 36).');

ro = [2.9 6.1 9.5 13.2 16.5 22.0];
ni = [48 22 14 10 8 6];
nt = [1036 709 535 421 352 278];
for ii = 1:6
    % seg_name = sprintf('seg_ni%d.nrrd', ni(ii));
    % v = nrrdread(seg_name);
    % v = v(:, :, end).';
    % if ii == 1
    %     v = double((v == 1));
    % else
    %     v = double((v == 2));
    % end

    t1_wm_subspace_std(ii) = mean_of_not_zero(v .* t1map_subspace_std(:, :, ii));
    t1_wm_maxgirf_std(ii) = mean_of_not_zero(v .* t1map_maxgirf_std(:, :, ii));
    t2_wm_subspace_std(ii) = mean_of_not_zero(v .* t2map_subspace_std(:, :, ii));
    t2_wm_maxgirf_std(ii) = mean_of_not_zero(v .* t2map_maxgirf_std(:, :, ii));

    % load(sprintf('./subspace_3d/pulseq_mrf_subspace_slice36_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii)), 't1map', 't2map');
    % t1_wm_subspace_std(ii) = std_of_not_zero(v .* t1map.');
    % t2_wm_subspace_std(ii) = std_of_not_zero(v .* t2map.');
    % figure;
    % imagesc(v .* t1map.', [0 1.5e3]);
    % axis image;
    % title(sprintf('T_{read} = %g ms', ro(ii)))
    % 
    % load(sprintf('./maxgirf_3d/pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii)), 't1map', 't2map');
    % t1_wm_maxgirf_std(ii) = std_of_not_zero(v .* t1map.');
    % t2_wm_maxgirf_std(ii) = std_of_not_zero(v .* t2map.');
end

f8_p2 = figure('Position', [380 300 780 530], 'Color', 'w');
yyaxis left;
plot(ro, t1_wm_subspace_std, '--', 'LineWidth', 4); hold on;
plot(ro, t1_wm_maxgirf_std, '-', 'LineWidth', 4);
plot(ro(5), t1_wm_maxgirf_std(5), 'Marker', '+', 'MarkerSize', 14, 'LineWidth', 4);
ylim([0 30]);
ylabel('T_1 Standard deviations [ms]', 'FontSize', 18);
yyaxis right;
plot(ro, t2_wm_subspace_std, '--', 'LineWidth', 4);
plot(ro, t2_wm_maxgirf_std, '-', 'LineWidth', 4);
plot(ro(5), t2_wm_maxgirf_std(5), 'Marker', '+', 'MarkerSize', 14, 'LineWidth', 4);
ylim([0 6]);
xlim([0 25]);

set(gca, 'FontSize', 18);
xlabel('T_{read} [ms]', 'FontSize', 18);
ylabel('T_2 Standard deviations [ms]', 'FontSize', 18);
% l = legend('gridding+subspace', 'MaxGIRF+subspace');

%% Save figure.
saveas(f8_p2, 'figure8_p2.fig');
saveas(f8_p2, 'figure8_p2.png');