clear; clc; close all;

indx = 3635; % T1 = 500 ms, T2 = 80 ms, WM.
indx2 = 275; % T1 = 100 ms, T2 = 100ms, ACR.

ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];

tr = [7 8.1 9.1 10.2 11.1 12.4 13.6 15.1 16.1 17.3 18.8 20.6 23 26.1];
ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];

%%
c = cool(14);
figure('Color', 'w');
hold on;
xlabel('Time [sec]');
xlim([0 7.5]);
set(gca, 'FontSize', 16);

for ii = 14:-1:1
    dict_name = sprintf('../../mrf/dictionaries/pulseq_readout_experiments/fisp_mrf_3d_ni%d_nt%d_rf1.mat', ni(ii), nt(ii));

    load(dict_name, 'dict', 'r');
    s = abs(dict(:, end, indx2)); % Signal evolution.
    mag_avg(15-ii) = mean(s);

    plot(((0:nt(ii)-1) * tr(15-ii) + 1.4) / 1000, s, 'Color', c(ii, :), 'LineWidth', 1.5);
end

%%
eff = sqrt(ro ./ tr);
eff = eff ./ eff(1);

figure('Color', 'w');
% colororder({'k', '#D95319'});

% figure('Color', 'w');
hold on;
xlabel('Readout, [ms]');
set(gca, 'FontSize', 16);

% yyaxis right;
plot(ro, eff, 'LineWidth', 0.8, 'Marker', 'x', 'Color', '#D95319', 'LineStyle', '--');
% ylabel('Relative SNR efficiency [unitless]');
% yyaxis left;
plot(ro, mag_avg / mag_avg(1), 'LineWidth', 1.5, 'Marker', '+', 'Color', '#0072BD');
% ylabel('Relative average magnitude [a.u.]');

% for ii = 14:-1:1
%     dict_name = sprintf('../mrf/dictionaries/pulseq_readout_experiments/fisp_mrf_3d_ni%d_nt%d_rf1.mat', ni(ii), nt(ii));
% 
%     load(dict_name, 'dict', 'r');
%     s = abs(dict(:, end, indx)); % Signal evolution.
%     mag_avg2(15-ii) = mean(s);
% end
plot(ro, mag_avg2 / mag_avg2(1), 'LineWidth', 1.5, 'Marker', '*', 'LineStyle', '-', 'Color', '#EDB120');

plot(ro, eff .* mag_avg / mag_avg(1), 'LineWidth', 1.5, 'Marker', '+', 'LineStyle', '--', 'Color', '#0072BD');
plot(ro, eff .* mag_avg2 / mag_avg2(1), 'LineWidth', 1.5, 'Marker', '*', 'LineStyle', '--', 'Color', '#EDB120');

ylabel('Relative values [uniteless]')
l = legend({'SNR efficiency', 'ACR signal', 'WM signal', 'ACR SNR', 'WM SNR'}, 'Location', 'northwest');
% legend('T_1=100ms, T_2=10ms', 'T_1=500ms, T_2=80ms');
