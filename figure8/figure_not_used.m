clear; clc; close all;

data_path = '/Users/zhibozhu/Documents/mrel_data/mrf/volunteer/readouts/pulseq_mrf_20240329';
addpath(data_path);

load(fullfile(data_path, 'mrf_vol1_analysis_slice36.mat'));
% v = nrrdread('seg_mprage_wm.nrrd');
% v = double(v(:, :, 36).');

ro = [2.9 6.1 9.5 13.2 16.5 22.0];
ni = [48 22 14 10 8 6];
nt = [1036 709 535 421 352 278];
for ii = 1:6
    seg_name = fullfile(data_path, sprintf('brain_ni%d.nrrd', ni(ii)));
    v = nrrdread(seg_name);
    v = v(:, :, end).';
    % if ii == 1
    %     v = double((v == 1));
    % else
    %     v = double((v == 2));
    % end
    v = (v == 2);

    load(fullfile(data_path, sprintf('subspace_3d/pulseq_mrf_subspace_slice36_B1__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), 't1map', 't2map');
    t1_brain_subspace = thre(not_zero(v .* t1map.'), 1e3);
    t2_brain_subspace = thre(not_zero(v .* t2map.'), 1e2);

    load(fullfile(data_path, sprintf('maxgirf_3d/pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii))), 't1map', 't2map');
    t1_brain_maxgirf = thre(not_zero(v .* t1map.'), 1e3);
    t2_brain_maxgirf = thre(not_zero(v .* t2map.'), 1e2);

    % figure;
    %     imagesc(thre(v .* t2map.', 1e2), [0 2e2]);
    % axis image;
    % title(sprintf('T_{read} = %g ms', ro(ii)));

    % figure;
    % histogram(t1_brain_subspace); hold on;
    % histogram(t1_brain_maxgirf); hold on;
    % title(sprintf('T_{read} = %g ms', ro(ii)));
end

f10 = figure('Position', [380 300 780 530], 'Color', 'w');
plot(ro, t1_brain_subspace, 'r--', 'LineWidth', 2); hold on;
plot(ro, t1_brain_maxgirf, 'r-', 'LineWidth', 2);
plot(ro, t2_brain_subspace, 'b--', 'LineWidth', 2);
plot(ro, t2_brain_maxgirf, 'b-', 'LineWidth', 2);
xlim([0 25]);

set(gca, 'FontSize', 14);
xlabel('T_{read} [ms]', 'FontSize', 18);
ylabel('Standard deviations [ms]', 'FontSize', 14, 'FontSize', 18);
legend('SD(T_1), before correction', 'SD(T_1), after correction', 'SD(T_2), before correction', 'SD(T_2), after correction');

