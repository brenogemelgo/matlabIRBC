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

        mon = char(tk);
        mon = replaceVars(mon);
        mon = strrep(mon, '*', ' * ');

        if isempty(symvar(ck))
            neg = isAlways(ck < 0);
            ak = simplify(abs(ck));

            if strcmp(mon, '1')
                term = string(castNum(ak));
            else

                if isequal(ak, sym(1))
                    term = string(mon);
                else
                    term = string(castNum(ak)) + " * " + string(mon);
                end

            end

        else
            ckStr = strtrim(char(ck));
            neg = startsWith(ckStr, '-');

            if neg
                ckStr = strtrim(ckStr(2:end));
            end

            ckStr = replaceVars(ckStr);
            ckStr = strrep(ckStr, '*', ' * ');

            if strcmp(mon, '1')
                term = string(ckStr);
            else

                if strcmp(ckStr, '1')
                    term = string(mon);
                else
                    term = "(" + string(ckStr) + ") * " + string(mon);
                end

            end

        end

        if first
            out = ternary(neg, "-" + term, term);
            first = false;
        else
            out = out + ternary(neg, " - " + term, " + " + term);
        end

    end

    s = char(out);
end

function r = ternary(cond, a, b)
    if cond, r = a; else, r = b; end
end
