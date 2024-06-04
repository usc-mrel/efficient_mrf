% ==== Show simulation of M(time,position,freq)
%
%	-> Try changing mode from 3 to 2...
%

clear;
close all;

b1 = [zeros(1,500) msinc(250,2) zeros(1,500)];
g = [zeros(1,375) -ones(1,125) ones(1,250) -ones(1,125) zeros(1,375)];
b1 = 1.57*b1/max(b1)/4;
x = [-4:.1:4];
T = .000004;
f = [-250:5:250];
t = [1:length(b1)]*T;
[mx,my,mz] = bloch(b1,g,t,1,.2,f,x,3); mxy=mx+i*my;
plot(t,[mx(:,41,35) my(:,41,35) mz(:,41,35)]);
xlabel('Time (s)');
ylabel('Magnetization');

