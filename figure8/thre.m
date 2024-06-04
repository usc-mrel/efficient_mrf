function out = thre(in, th)
    [nx, ny] = size(in);
    indx = abs(in) >= th;
    out = in;
    out(indx) = 0;
end