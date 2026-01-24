function s = polynom(P, vars)
    P = expand(P);

    if isequal(P, sym(0))
        s = 'static_cast<scalar_t>(0)';
        return;
    end

    [c, t] = coeffs(P, vars);

    out = "";
    first = true;

    for k = 1:numel(c)
        ck = simplify(c(k));
        tk = t(k);

        if isequal(ck, sym(0))
            continue;
        end

        neg = isAlways(ck < 0);
        ak = simplify(abs(ck));

        mon = char(tk);
        mon = replaceVars(mon);
        mon = strrep(mon, '*', ' * ');

        if strcmp(mon, '1')
            term = string(castNum(ak));
        else

            if isequal(ak, sym(1))
                term = string(mon);
            else
                term = string(castNum(ak)) + " * " + string(mon);
            end

        end

        if first

            if neg
                out = "-" + term;
            else
                out = term;
            end

            first = false;
        else

            if neg
                out = out + " - " + term;
            else
                out = out + " + " + term;
            end

        end

    end

    s = char(out);
end
