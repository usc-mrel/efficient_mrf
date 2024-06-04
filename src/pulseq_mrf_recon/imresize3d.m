function img2 = imresize3d(img,s,method)
%   Interpolates images in 3D

%   Check inputs
if nargin < 3 || isempty(method)
    method = 'linear';
end
if nargin < 2
    s = 2;
end
if nargin < 1
    error('Requires at least one input');
end

%   Get image size
[np nv ns] = size(img);
if ns ~= size(img,3);
    error('Image must be of size NP x NV x NS');
end

%   Check desired interpolation
if numel(s) == 1
    s = s.*[np nv ns];
elseif numel(s) ~= 3
    error('Wrong interpolation size');
end

%   Create interpolation points
[x y z] = ndgrid(1:np,1:nv,1:ns);
x = single(x);y = single(y);z = single(z);
[X Y Z] = ndgrid(linspace(1,np,s(1)),...
                 linspace(1,nv,s(2)),...
                 linspace(1,ns,s(3)));
X = single(X);Y = single(Y);Z = single(Z);

%   Interpolate
if strcmp(method,'linear')
    img2 = interpn(x,y,z,img,X,Y,Z,'linear');
elseif strcmp(method,'spline')
    img2 = interpn(x,y,z,img,X,Y,Z,'spline');
elseif strcmp(method,'cubic')
    img2 = interpn(x,y,z,img,X,Y,Z,'cubic');
else
    img2 = interpn(x,y,z,img,X,Y,Z,'nearest');
end
