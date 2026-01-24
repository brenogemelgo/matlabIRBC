function tf = exprEqual(a, b)
    d = simplify(a - b, 'Steps', 200);

    if isequal(d, sym(0))
        tf = true;
        return;
    end

    if isempty(symvar(d))
        tf = false;
        return;
    end

    vars = [sym('rhoI'), sym('mxyI'), sym('mxzI'), sym('myzI'), sym('omega')];
    omegaVals = [sym(2) / 3, sym(3) / 5, sym(5) / 7, sym(7) / 11, sym(11) / 13];
    ints = sym([2 3 5 7 11 13 17 19 23 29]);

    hits = 0;

    for t = 1:10
        vals = [ints(t), ints(mod(t, 10) + 1), ints(mod(t + 1, 10) + 1), ints(mod(t + 2, 10) + 1), omegaVals(mod(t - 1, numel(omegaVals)) + 1)];
        dt = simplify(subs(d, vars, vals));

        [~, den] = numden(dt);

        if isequal(simplify(den), sym(0))
            continue;
        end

        hits = hits + 1;

        if ~isequal(dt, sym(0))
            tf = false;
            return;
        end

        if hits >= 4
            tf = true;
            return;
        end

    end

    tf = false;
end
