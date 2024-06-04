function visualizeMagn(RF, B0, M, accel, filename)

% visualizeMagn(RF, B0, M)
%
% function visualizes the magnetization for a given RF B0 and magnetization
% 
% Input: 
%       RF is [N,1] representing B1 in [G] (real is y axis, imag is x axis)
%       B0 is [N,1] The effecive B0 field due to off-resonance or gradient
%       M  is [N,3] representing Mx, My and Mz as a function of time
%       accel - factor of acceleration of plot 
% (c) 2011 Michael Lustig

if nargin < 5
    filename = [];
end


ppt_font = 'Comic Sans MS';
ppt_size = 32;
ppt_lw = 4;
ppt_col = 'w';

fig = figure;
hold on;
% Draw axes
%
sc = 1/0.85;
x0 = -sc * eye(3);
x1 = sc * eye(3);
for ix = 1:3, 
    plot3([x0(ix,1) x1(ix,1)], [x0(ix,2) x1(ix,2)], [x0(ix,3) x1(ix,3)],'w', ...
          'LineWidth', ppt_lw/2);
end;

maxB = max([abs(B0(:));abs(RF(:))]);
RF = RF/maxB;
B0 = B0/maxB;


axis off
view([4 2 2]);
lighting phong;
camlight right;
zoom(2.2);

h1 =[];
h2 = [];
h3 = [];
h4 = [];

if isempty(filename) == 0
    aviobj = avifile(filename)
end

for n = 1:accel:length(RF)
    delete(h1);
    h1 = arrow3D([0 0 0],[real(RF(n)) imag(RF(n)) 0],'y',0.8,0.03);
    delete(h2);
    h2 = arrow3D([0 0 0],[0 0 B0(n)],'c',0.8,0.03);
    delete(h3);
    h3 = arrow3D([0 0 0],[M(n,1) M(n,2) M(n,3)],'r',0.8,0.03);
    delete(h4);
    h4 = plot3(M(1:n,1),M(1:n,2),M(1:n,3),'k--','LineWidth', ppt_lw/2);
    drawnow;
    if isempty(filename) == 0
        F = getframe(fig);
        aviobj = addframe(aviobj,F);
    end
    
end

n = length(RF);
delete(h1);
    h1 = arrow3D([0 0 0],[real(RF(n)) imag(RF(n)) 0],'y',0.8,0.03);
    delete(h2);
    h2 = arrow3D([0 0 0],[0 0 B0(n)],'c',0.8,0.03);
    delete(h3);
    h3 = arrow3D([0 0 0],[M(n,1) M(n,2) M(n,3)],'r',0.8,0.03);
    delete(h4);
    h4 = plot3(M(1:n,1),M(1:n,2),M(1:n,3),'k--','LineWidth', ppt_lw/2);
    drawnow;
    if isempty(filename)==0
        F = getframe(fig);
        aviobj = addframe(aviobj,F);
        aviobj = close(aviobj);        
    end
end