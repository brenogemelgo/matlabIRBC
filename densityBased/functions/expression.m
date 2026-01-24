function s = expression(expr)
    expr = simplify(expr);

    if isequal(expr, sym(0))
        s = 'static_cast<scalar_t>(0)';
        return;
    end

    vars = [sym('rhoI'), sym('mxyI'), sym('mxzI'), sym('myzI'), sym('omega')];

    [N, D] = numden(expr);
    N = expand(N);
    D = expand(D);

    L = sym(1);
    [cN, ~] = coeffs(N, vars);
    [cD, ~] = coeffs(D, vars);
    cAll = [cN, cD];

    for k = 1:numel(cAll)
        [~, den] = numden(cAll(k));
        L = lcm(L, den);
    end

    Nint = expand(L * N);
    Dint = expand(L * D);

    numStr = polynom(Nint, vars);
    denStr = polynom(Dint, vars);

    if isequal(Dint, sym(1))
        s = numStr;
    else
        s = ['(' numStr ') / (' denStr ')'];
    end

end
