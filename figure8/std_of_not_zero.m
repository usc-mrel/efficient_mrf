function out = std_of_not_zero(in)
    in = in(abs(in) ~= 0);
    out = std(in);
end