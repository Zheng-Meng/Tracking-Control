function [rmse] = func_rmse(a, b, time_start, time_end)
[c,d] = size(a);
len = max(c, d);

if size(a, 2) ~= len
    a = a';
end

if size(b, 2) ~= len
    b = b';
end

rmse = sqrt(mean( sum( (a(:, time_start:time_end)-b(:, time_start:time_end)).^2, 1 )));

end