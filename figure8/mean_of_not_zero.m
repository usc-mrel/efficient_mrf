function out = mean_of_not_zero(in)
    in = in(abs(in) ~= 0);
    out = mean(in);
end