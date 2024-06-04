%
%	function [k, g, s, m1, m2, t, v] = plotgradinfo(g,T,k0)
%
%	Function makes a nice plot given a 2D gradient waveform.
%
%	INPUT:
%		g	gradient (G/cm) = gx + i*gy, or Nx2, Nx3 or Nx4.
%		T	sample period (s)
%		k0	initial condition for k-space.
%
%	OUTPUT:
%		k	k-space trajectory (cm^(-1))
%		g	gradient (G/cm)
%		s	slew rate trajectory (G/cm/s)
%		m1	first moment trajectory (s/cm)
%		m2	second moment trajectory (s^2/cm)
%
%	B.Hargreaves, Aug 2002.
%

% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: plotgradinfo.m,v $
%	Revision 1.7  2003/09/16 16:45:03  brian
%	Added ability to plot acquisition interval.
%	
%	Revision 1.6  2003/09/03 18:42:32  brian
%	Added a comment.  Voltage to 2nd moment.
%	
%	Revision 1.5  2002/10/14 15:59:36  brian
%	Added calculation/plot of voltage.
%	
%	Revision 1.4  2002/09/05 18:35:44  bah
%	Separted into calculation and plot parts.
%	
%	Revision 1.3  2002/09/05 17:45:31  bah
%	Plots and returns first and second moments now too.
%	
%	Revision 1.2  2002/09/04 01:20:17  bah
%	Fixed scaling of slew, and other bugs.
%	
%	Revision 1.1  2002/09/03 23:08:20  bah
%	Added to CVS.
%	
%
% ===========================================================


function [k, g, s, m1, m2, t, v] = plotgradinfo(g, T, k0)


if (nargin < 2)
    T = .000004; % [sec]
end
if (nargin < 3)
    k0 = 0;
end
gamma = 4258; % [Hz/G]

sg = size(g);
if ((sg(2)==1) && ~isreal(g))
	g = [real(g) imag(g)];
end

[k, g, s, m1, m2, t, v] = calcgradinfo(g, T, k0);

sg = size(g);
ng = sg(2);
if (ng > 3)
	acq = g(:,4);
	g = g(:,1:3);
	k = k(:,1:3);
	s = s(:,1:3);
	m1 = m1(:,1:3);
	m2 = m2(:,1:3);
	v = v(:,1:3);
	sg = size(g);
	ng = sg(2);
else
	acq = [];
end

%pls={'g--','r--','b--','k:'};
pls = {[31 119 180]/255, [255 127 14]/255, [44 160 44]/255, [214 39 40]/255};
%pls{ng+1}='k:';

M=2;
N=3;

clf;
set(gcf, 'Color', 'w');
if (ng >= 2)
    subplot(M,N,1);
    if (max(size(acq)) > 0)
        mxk = sqrt(max(acq'.*sum((k(:,1:2))'.*(k(:,1:2))')));
        xsq = cos(2 * pi * (0:30).' / 30) * mxk;
        ysq = sin(2 * pi * (0:30).' / 30) * mxk;
        h = patch(xsq, ysq, .8 * [1 1 1]);
        set(h, 'EdgeColor', .8 * [1 1 1]);
        hold on;
    end
    plot(k(:,1), k(:,2), 'Color', pls{1}, 'LineWidth', 1);
    xlabel('k_x [cycle/cm]');
    ylabel('k_y [cycle/cm]');
    title('K-space Trajectory');
    grid on;
    a = axis;
    axis(max(abs(a))*[-1 1 -1 1]);
end

subplot(M, N, 2);
doplot(t * 1e3, k, 'k', pls, acq);
xlabel('Time [msec]');
ylabel('K-space [cycle/cm]');
title('K-space vs time');

subplot(M, N, 3);
doplot(t * 1e3, g, 'g', pls, acq);
xlabel('Time [msec]');
ylabel('Gradient [G/cm]');
title('Gradient vs time');

subplot(M, N, 4);
doplot(t * 1e3, s, 'slew', pls, acq);
xlabel('Time [msec]');
ylabel('Slew Rate [G/cm/s]');
title('Slew Rates vs time');

if (M * N >= 5)
    subplot(M, N, 5);
    doplot(t * 1e3, m1, 'm1', pls, acq);
    xlabel('Time [msec]');
    ylabel('First Moment [s/cm]');
    title('First Moment vs time');
end

if (0==0) % Plot second moment
    if (M *N >= 6)
        subplot(M, N, 6);
        doplot(t * 1e3, m2,' m2', pls, acq);
        xlabel('Time [msec]');
        ylabel('Second Moment [s^2/cm]');
        title('Second Moment vs time');
    end
else % Plot coil voltage.
    if (M * N >= 6)
        subplot(M, N, 6);
        doplot(t * 1e3, v, 'voltage', pls, acq);
        xlabel('Time [msec]');
        ylabel('Voltage [V]');
        title('Coil Voltage vs time');
    end
end


function doplot(x, y, labs, col, acq)

axl = {'x', 'y', 'z'};
ss = size(y);
if (min(size(acq))>0)
    mina = min(find(acq > 0));
    maxa = max(find(acq > 0));
    xsq = [x(mina) x(maxa) x(maxa) x(mina)];
    ysq = 1.1*sqrt(max(abs(sum(y'.*y'))))*[-1 -1 1 1];
    patchc = .8*[1 1 1];
end
hold off;
if (min(size(acq)) > 0)
    h = patch(xsq, ysq, patchc);
    set(h, 'EdgeColor', patchc);
    hold on;
end
for q = 1:ss(2)
    plot(x, y(:,q), 'Color', col{q}, 'LineWidth', 1);
    hold on;
    leg{q} = sprintf('%s_%s', labs, axl{q});
end
if (ss(2) > 1)
    ab = sqrt(sum(y.'.*y.')).';
    plot(x, ab, 'Color', col{ss(2)+1}, 'LineWidth', 1);
    leg{ss(2) + 1} = sprintf('|%s|', labs);
end
hold off;
legend(leg, 'Location', 'northwest');
grid on;
