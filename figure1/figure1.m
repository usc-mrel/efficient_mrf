clear; clc; close all

%%% Tight subplot path.
addpath(genpath('../src/utils/tight_subplot'));

ni = [6 7 8 9 10 11 12 14 16 19 22 27 34 48];
nt = [278 317 352 382 421 451 481 535 585 654 709 798 892 1036];

%%% Load spatial coordinates.
load('concomitant_field_ni30.mat', 'p');

tr = [7 8.1 9.1 10.2 11.1 12.4 13.6 15.1 16.1 17.3 18.8 20.6 23 26.1];
ro = [2.9 4.0 5.0 6.1 7.0 8.3 9.5 11.1 12.1 13.2 14.7 16.5 18.9 22.0];

%% Figure 1, spiral arms and concomitant field induced phases.
gamma = 4257.59 * (1e4 * 2 * pi); % [rad/sec/T]
Nl = 19;
B0 = 0.55;
dz = [0 0.05 0.1 0.15];

clear phase;
f1 = figure('Position', [260 500 910 590], 'Color', 'w', 'InvertHardcopy', 'off');
text_color = [0.85 0.85 0.85];
[ha, pos] = tight_subplot(2, 3, [0.1 0.1], [0.1 0.08], [0.08 0.05]);
ha = (reshape(ha, [3 2])).';

%%% RO = 2.9 ms.
load(sprintf('traj_ni%d_nt%d.mat', ni(14), nt(14)));
krmax = 2 * pi / (0.3 / 256) /2;
k_rcs = kx_GIRF(:, 1) * 2 * krmax;
k_rcs(:, :, 2) = -ky_GIRF(:, 1) * 2 * krmax;
k_rcs(:, :, 3) = 0;
k_rcs = permute(k_rcs, [3 1 2]);
R_rcs2gcs = [0 1 0; 1 0 0; 0 0 1];
k_gcs = R_rcs2gcs * k_rcs;
g_gcs = diff(k_gcs, [], 2) / (4257.59e4 * 2 * pi * 1e-3 * 2.5e-6); % [mT/m]
R_gcs2dcs = [0 1 0; -1 0 0; 0 0 -1];
g_dcs = R_gcs2dcs * g_gcs * 1e-1;
k = calculate_concomitant_field_coefficients(permute(squeeze(g_dcs(1, :)), [2 1]), permute(squeeze(g_dcs(2, :)), [2 1]), permute(squeeze(g_dcs(3, :)), [2 1]), Nl, B0, gamma, dt);

for ii = 1:4
    p(:, 3) = dz(ii);
    p(:, 6) = p(:, 3).^2;
    temp = k(:, 4:end) * p(:, 4:end).';
    phase(:, ii) = mean(temp, 2) / (2 * pi);
end

axes(ha(1, 1));
plot(kx_GIRF(:, 1), ky_GIRF(:, 1), 'LineWidth', 2);
xlabel('k_x', 'FontSize', 16);
ylabel('k_y', 'FontSize', 16);
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
axis square;
title(sprintf('RO = %g ms', ro(1)), 'FontSize', 24);

axes(ha(2, 1));
dz = [0 0.05 0.1 0.15];
for ii = 1:4
    p(:, 3) = dz(ii);
    p(:, 6) = p(:, 3).^2;
    phase = k(:, 4:end) * p(:, 4:end).';
    phase = mean(phase, 2) / (2 * pi);
    plot(dt * (0:size(phase, 1)) * 1e3, [0; phase], 'LineWidth', 2);
    hold on;
end
xlabel('Time [ms]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Phase [cycles]', 'FontSize', 16, 'FontWeight', 'bold');
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlim([0 17]);
ylim([0 8]);
legend('|z|=0mm', '|z|=50mm', '|z|=100mm', '|z|=150mm', 'Location', 'northwest', 'FontSize', 14);
clear phase;

%%% RO = 9.5 ms.
load(sprintf('traj_ni%d_nt%d.mat', ni(8), nt(8)));
krmax = 2 * pi / (0.3 / 256) /2;
k_rcs = kx_GIRF(:, 1) * 2 * krmax;
k_rcs(:, :, 2) = -ky_GIRF(:, 1) * 2 * krmax;
k_rcs(:, :, 3) = 0;
k_rcs = permute(k_rcs, [3 1 2]);
R_rcs2gcs = [0 1 0; 1 0 0; 0 0 1];
k_gcs = R_rcs2gcs * k_rcs;
g_gcs = diff(k_gcs, [], 2) / (4257.59e4 * 2 * pi * 1e-3 * 2.5e-6); % [mT/m]
R_gcs2dcs = [0 1 0; -1 0 0; 0 0 -1];
g_dcs = R_gcs2dcs * g_gcs * 1e-1;
k = calculate_concomitant_field_coefficients(permute(squeeze(g_dcs(1, :)), [2 1]), permute(squeeze(g_dcs(2, :)), [2 1]), permute(squeeze(g_dcs(3, :)), [2 1]), Nl, B0, gamma, dt);

for ii = 1:4
    p(:, 3) = dz(ii);
    p(:, 6) = p(:, 3).^2;
    temp = k(:, 4:end) * p(:, 4:end).';
    phase(:, ii) = mean(temp, 2) / (2 * pi);
end

axes(ha(1, 2));
plot(kx_GIRF(:, 1), ky_GIRF(:, 1), 'LineWidth', 2);
% text(-3, 2.75, sprintf('RO=%gms', ro(15 - 14)), 'FontSize', 12, 'HorizontalAlignment', 'left');
% text(-3, 2.4, sprintf('TR=%gms', tr(15 - 14)), 'FontSize', 12, 'HorizontalAlignment', 'left');
xlabel('k_x', 'FontSize', 16);
ylabel('k_y', 'FontSize', 16);
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
axis square;
title(sprintf('RO = %g ms', ro(7)), 'FontSize', 24);

axes(ha(2, 2));
dz = [0 0.05 0.1 0.15];
for ii = 1:4
    p(:, 3) = dz(ii);
    p(:, 6) = p(:, 3).^2;
    phase = k(:, 4:end) * p(:, 4:end).';
    phase = mean(phase, 2) / (2 * pi);
    plot(dt * (0:size(phase, 1)) * 1e3, [0; phase], 'LineWidth', 2);
    hold on;
end
xlabel('Time [ms]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Phase [cycles]', 'FontSize', 16, 'FontWeight', 'bold');
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlim([0 17]);
ylim([0 8]);
legend('|z|=0mm', '|z|=50mm', '|z|=100mm', '|z|=150mm', 'Location', 'northwest', 'FontSize', 14);
clear phase;

%%% RO = 16.5 ms.
load(sprintf('traj_ni%d_nt%d.mat', ni(3), nt(3)));
krmax = 2 * pi / (0.3 / 256) /2;
k_rcs = kx_GIRF(:, 1) * 2 * krmax;
k_rcs(:, :, 2) = -ky_GIRF(:, 1) * 2 * krmax;
k_rcs(:, :, 3) = 0;
k_rcs = permute(k_rcs, [3 1 2]);
R_rcs2gcs = [0 1 0; 1 0 0; 0 0 1];
k_gcs = R_rcs2gcs * k_rcs;
g_gcs = diff(k_gcs, [], 2) / (4257.59e4 * 2 * pi * 1e-3 * 2.5e-6); % [mT/m]
R_gcs2dcs = [0 1 0; -1 0 0; 0 0 -1];
g_dcs = R_gcs2dcs * g_gcs * 1e-1;
k = calculate_concomitant_field_coefficients(permute(squeeze(g_dcs(1, :)), [2 1]), permute(squeeze(g_dcs(2, :)), [2 1]), permute(squeeze(g_dcs(3, :)), [2 1]), Nl, B0, gamma, dt);

for ii = 1:4
    p(:, 3) = dz(ii);
    p(:, 6) = p(:, 3).^2;
    temp = k(:, 4:end) * p(:, 4:end).';
    phase(:, ii) = mean(temp, 2) / (2 * pi);
end

axes(ha(1, 3));
plot(kx_GIRF(:, 1), ky_GIRF(:, 1), 'LineWidth', 2);
% text(-3, 2.75, sprintf('RO=%gms', ro(15 - 14)), 'FontSize', 12, 'HorizontalAlignment', 'left');
% text(-3, 2.4, sprintf('TR=%gms', tr(15 - 14)), 'FontSize', 12, 'HorizontalAlignment', 'left');
xlabel('k_x', 'FontSize', 16);
ylabel('k_y', 'FontSize', 16);
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
axis square;
title(sprintf('RO = %g ms', ro(12)), 'FontSize', 24);

axes(ha(2, 3));
dz = [0 0.05 0.1 0.15];
for ii = 1:4
    p(:, 3) = dz(ii);
    p(:, 6) = p(:, 3).^2;
    phase = k(:, 4:end) * p(:, 4:end).';
    phase = mean(phase, 2) / (2 * pi);
    plot(dt * (0:size(phase, 1)) * 1e3, [0; phase], 'LineWidth', 2);
    hold on;
end
xlabel('Time [ms]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Phase [cycles]', 'FontSize', 16, 'FontWeight', 'bold');
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlim([0 17]);
ylim([0 8]);
legend('|z|=0mm', '|z|=50mm', '|z|=100mm', '|z|=150mm', 'Location', 'northwest', 'FontSize', 14);
clear phase;

%% Save figures.
saveas(f1, 'figure1.fig');
saveas(f1, 'figure1.png');


