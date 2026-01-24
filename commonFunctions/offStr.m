function s = offStr(d)

    if d == 0
        s = '';
    elseif d > 0
        s = sprintf(' + %d', d);
    else
        s = sprintf(' - %d', -d);
    end

end