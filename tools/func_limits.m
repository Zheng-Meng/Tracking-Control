function [value] = func_limits(value, limit_type)
%% limitations

qdt_lim = [-0.05, 0.05];
q2dt_lim = [-0.3, 0.3];

% limit type:
%   1 : dq/dt
%   2 : d2q/dt2

if limit_type == 1
    if value(1) > qdt_lim(2)
        value(1) = qdt_lim(2);
    end
    if value(1) < qdt_lim(1)
        value(1) = qdt_lim(1);
    end
    if value(2) > qdt_lim(2)
        value(2) = qdt_lim(2);
    end
    if value(2) < qdt_lim(1)
        value(2) = qdt_lim(1);
    end
elseif limit_type == 2
    if value(1) > q2dt_lim(2)
        value(1) = q2dt_lim(2);
    end
    if value(1) < q2dt_lim(1)
        value(1) = q2dt_lim(1);
    end
    if value(2) > q2dt_lim(2)
        value(2) = q2dt_lim(2);
    end
    if value(2) < q2dt_lim(1)
        value(2) = q2dt_lim(1);
    end
end

end

