% demo_Turley_2013_MRM.m
% Written by Namgyun Lee
% Email: namgyunl@kbsi.re.kr, ggang56@gmail.com (preferred)
% Started: 12/27/2016, Last modified: 12/27/2016
% Updated by Nam Gyun Lee, namgyunl@usc.edu, 12/18/2022

%% Clean slate
close all; clear all; clc;

%%
% Units are in Hz, sec, gauss, and cm!

%% Define input parameters
gamma    = 4257.746778; % typically 4257 [Hz/G]
fov      = 24.0;        % field of view [cm]
res      = 0.1;         % resolution [cm]
slewmax  = 130e2;       % max slew rate [G/cm/sec]
gsysmax  = 4;           % total max gradient amplitude [G/cm]
rast     = 4e-6;        % "raster time" between samples in sec, i.e., 0.000004 gives 4 usec gradient update
arms     = 25;          % number of spiral interleaves
greadmax = min(gsysmax, 2 /(gamma * fov * rast)); % max gradient amplitude (system or bw-limited) [G/cm]

%--------------------------------------------------------------------------
% The next 4 variables are for variable density spirals
% They create a transition in the radial spacing as the kspace radius goes
% from 0 to 1, i.e.,
%     0 < kr < us_0, spacing = nyquist distance
%  us_0 < kr < us_1, spacing increases to us_r (affected by ustype)
%  us_1 < kr < 1   , spacing = us_r
% SET us_0 = 1 or us_r = 1 to avoid variable sampling
%
% rate of change in undersampling (see code below)
% ustype: 0 = linearly changing undersampling
%         1 = quadratically changing undersampling
%         2 = Hanning-type chnage in undersampling
%--------------------------------------------------------------------------
us_0 = 1;
us_1 = 1;
us_r = 1;
ustype = 0;

%--------------------------------------------------------------------------
% Some more variables determining the type of waveform
% For gtype: 0 = calculate through readout
%            1 = include grad rampdown
%            2 = include rewinder to end at k=0
%--------------------------------------------------------------------------
gtype  = 2;
sptype = 0; % 0 = Archimedean, 1 = Fermat

%% Calculate a base k-space trajectory with spiralgen.c
spparams = [gamma;
            fov;
            res;
            slewmax;
            greadmax;
            gsysmax;
            rast;
            arms;
            us_0;
            us_1;
            us_r;
            ustype;
            gtype;
            sptype];

%--------------------------------------------------------------------------
% This function takes parameters passed in spparams array and
% returns a single spiral arm calculated numerically.
% The k-space trajectory for the arm is in karray
% The corresponding gradient waveform is in garray
% The parameter garrlen are the lengths of karray and garray
% The parameter karrlen indicates when sampling stops
%
% FINAL ARRAYS are {kx(0), ky(0), kx(1), ky(1), ...} and
%                  {gx(0), gy(0), gx(1), gy(1), ...}
%--------------------------------------------------------------------------
[karray, garray, karrlen, garrlen] = spiralgen_jgp_11apr_mex(spparams);

% Gradient coordinate system: PRS = [PE,RO,SL] = [gu,gv,gw]
kv = karray(1:2:garrlen*2); % RO [cycle/cm]
ku = karray(2:2:garrlen*2); % PE [cycle/cm]
gv = garray(1:2:garrlen*2) * 1e1; % RO [G/cm] * [T/1e4G] * [1e3mT/T] * [1e2cm/m] => *1e1 [mT/m]
gu = garray(2:2:garrlen*2) * 1e1; % PR [G/cm] * [T/1e4G] * [1e3mT/T] * [1e2cm/m] => *1e1 [mT/m]

T_read = karrlen * rast;
T_all = garrlen * rast;

%% Calculate slew rate [mT/m/msec]
sv = cat(1, diff(gv,1), 0) / rast * 1e-3; % [mT/m] / [sec] * [sec/1e3msec] => *1e-3 [mT/m/msec]
su = cat(1, diff(gu,1), 0) / rast * 1e-3; % [mT/m] / [sec] * [sec/1e3msec] => *1e-3 [mT/m/msec]

%% Calculate a time axis [sec]
t = (0:garrlen-1).' * rast; % [sec]

%% Display k-space trajectories, gradient waveforms, and slew rates
FontSize = 14;

kvmax = max(abs(kv));
kumax = max(abs(ku));
kmax = max(kvmax, kumax);

figure('Color', 'w', 'Position', [1 245 1382 745]);
subplot(3,2,1);
hold on;
grid on; grid minor;
plot(t * 1e3, kv, 'LineWidth', 1);
plot(t * 1e3, ku, 'LineWidth', 1);
xlim([0 max(t)] * 1e3);
ylim([-kmax kmax] * 1.1);
set(gca, 'Box', 'On', 'FontSize', FontSize-2);
xlabel('Time [msec]', 'FontSize', FontSize-2);
ylabel('k-space [cycle/cm]', 'FontSize', FontSize-2);
legend('kv (RO)', 'ku (PE)', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', FontSize-4);

subplot(3,2,3);
hold on;
grid on; grid minor;
plot(t * 1e3, gv, 'LineWidth', 1);
plot(t * 1e3, gu, 'LineWidth', 1);
plot(t * 1e3, sqrt(gv.^2 + gu.^2));
xlim([0 max(t)] * 1e3);
set(gca, 'Box', 'On', 'FontSize', FontSize-2);
xlabel('Time [msec]', 'FontSize', FontSize-2);
ylabel('Gradient [mT/m]', 'FontSize', FontSize-2);
legend('Gv (RO)', 'Gu (PE)', '|G|', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', FontSize-4);

subplot(3,2,5);
hold on;
grid on; grid minor;
plot(t * 1e3, sv, 'LineWidth', 1);
plot(t * 1e3, su, 'LineWidth', 1);
plot(t * 1e3, sqrt(sv.^2 + su.^2));
xlim([0 max(t)] * 1e3);
set(gca, 'Box', 'On', 'FontSize', FontSize-2);
xlabel('Time [msec]', 'FontSize', FontSize-2);
ylabel('Slew rate [mT/m/msec]', 'FontSize', FontSize-2);
%ylim([-150 150]);
legend('Sv (RO)', 'Su (PE)', '|S|', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', FontSize-4);

subplot(3,2,[2 4 6]);
hold on;
grid on; grid minor;
plot(kv, ku, 'LineWidth', 1);
axis image;
set(gca, 'Box', 'On', 'FontSize', FontSize-4);
xlabel('kv (RO) [cycle/cm]', 'FontSize', FontSize-2);
ylabel('ku (PE) [cycle/cm]', 'FontSize', FontSize-2);
title({'Base k-space trajectory (readout + rewinder)', sprintf('Tread/Tall = %5.3f/%5.3f msec', T_read*1e3, T_all*1e3)});
xlim([-kmax kmax] * 1.1);
ylim([-kmax kmax] * 1.1);
