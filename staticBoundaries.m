clc; clearvars; close all

stencil = "D3Q27"; % "D3Q19" or "D3Q27"
phase_field = true;

[Q, w, cx, cy, cz, bitLists, bcs] = defineStencil(stencil);

% constants
as = sqrt(sym(3));
as2 = as ^ 2;
as4 = as ^ 4;
cxS = sym(cx(:));
cyS = sym(cy(:));
czS = sym(cz(:));
wS = sym(w(:));

% hermite tensors
Hxx = cxS .^ 2 - as ^ (-2);
Hxy = cxS .* cyS;
Hxz = cxS .* czS;
Hyy = cyS .^ 2 - as ^ (-2);
Hyz = cyS .* czS;
Hzz = czS .^ 2 - as ^ (-2);

% symbols
pI = sym('pI', 'real');
mxyI = sym('mxyI', 'real');
mxzI = sym('mxzI', 'real');
myzI = sym('myzI', 'real');
omega = sym('omega', 'real');

% cases
corner = { ...
              struct('name', "WEST_SOUTH_BACK", 'key', "SOUTH_WEST_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "WEST_SOUTH_FRONT", 'key', "SOUTH_WEST_FRONT", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_SOUTH_BACK", 'key', "SOUTH_EAST_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_SOUTH_FRONT", 'key', "SOUTH_EAST_FRONT", 'solve', string.empty(1, 0)), ...
              struct('name', "WEST_NORTH_BACK", 'key', "NORTH_WEST_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "WEST_NORTH_FRONT", 'key', "NORTH_WEST_FRONT", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_NORTH_BACK", 'key', "NORTH_EAST_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_NORTH_FRONT", 'key', "NORTH_EAST_FRONT", 'solve', string.empty(1, 0)) ...
          };

edge = { ...
            struct('name', "WEST_SOUTH", 'key', "SOUTH_WEST", 'solve', ["mxy"]), ...
            struct('name', "EAST_SOUTH", 'key', "SOUTH_EAST", 'solve', ["mxy"]), ...
            struct('name', "WEST_NORTH", 'key', "NORTH_WEST", 'solve', ["mxy"]), ...
            struct('name', "EAST_NORTH", 'key', "NORTH_EAST", 'solve', ["mxy"]), ...
            struct('name', "WEST_BACK", 'key', "WEST_BACK", 'solve', ["mxz"]), ...
            struct('name', "WEST_FRONT", 'key', "WEST_FRONT", 'solve', ["mxz"]), ...
            struct('name', "EAST_BACK", 'key', "EAST_BACK", 'solve', ["mxz"]), ...
            struct('name', "EAST_FRONT", 'key', "EAST_FRONT", 'solve', ["mxz"]), ...
            struct('name', "SOUTH_BACK", 'key', "SOUTH_BACK", 'solve', ["myz"]), ...
            struct('name', "SOUTH_FRONT", 'key', "SOUTH_FRONT", 'solve', ["myz"]), ...
            struct('name', "NORTH_BACK", 'key', "NORTH_BACK", 'solve', ["myz"]), ...
            struct('name', "NORTH_FRONT", 'key', "NORTH_FRONT", 'solve', ["myz"]) ...
        };

face = { ...
            struct('name', "WEST", 'key', "WEST", 'solve', ["mxy", "mxz"]), ...
            struct('name', "EAST", 'key', "EAST", 'solve', ["mxy", "mxz"]), ...
            struct('name', "SOUTH", 'key', "SOUTH", 'solve', ["mxy", "myz"]), ...
            struct('name', "NORTH", 'key', "NORTH", 'solve', ["mxy", "myz"]), ...
            struct('name', "BACK", 'key', "BACK", 'solve', ["mxz", "myz"]), ...
            struct('name', "FRONT", 'key', "FRONT", 'solve', ["mxz", "myz"]) ...
        };

allCases = [corner(:); edge(:); face(:)];
allCases = allCases.';

% main loop
for c = 1:numel(allCases)
    name = allCases{c}.name;
    key = allCases{c}.key;

    solveShears = string(allCases{c}.solve);
    solveShears = solveShears(:).';

    Os = getCoordinates(bcs.(key), bitLists);
    Is = getOppositeCoordinates(Os, cx, cy, cz);

    % incoming moments
    % fprintf('const scalar_t p_I = %s;\n', cpp_sum_pop(Is));
    % if any(solveShears == "mxy")
    %     fprintf('const scalar_t mxy_I = %s;\n', cpp_sum_pop_weighted(Is, Hxy(Is + 1)));
    % end
    % if any(solveShears == "mxz")
    %     fprintf('const scalar_t mxz_I = %s;\n', cpp_sum_pop_weighted(Is, Hxz(Is + 1)));
    % end
    % if any(solveShears == "myz")
    %     fprintf('const scalar_t myz_I = %s;\n', cpp_sum_pop_weighted(Is, Hyz(Is + 1)));
    % end

    % solve boundary moments
    sol = solve_static_case(Q, wS, Hxx, Hxy, Hxz, Hyy, Hyz, Hzz, Os, Is, as2, as4, omega, pI, mxyI, mxzI, myzI, solveShears);

    % print solved boundary moments
    % fprintf('const scalar_t p = %s;\n', cpp_expr(sol.p));
    % if isfield(sol, 'mxy')
    %     fprintf('const scalar_t mxy = %s;\n', cpp_expr(sol.mxy));
    % end
    % if isfield(sol, 'mxz')
    %     fprintf('const scalar_t mxz = %s;\n', cpp_expr(sol.mxz));
    % end
    % if isfield(sol, 'myz')
    %     fprintf('const scalar_t myz = %s;\n', cpp_expr(sol.myz));
    % end

    fprintf('\ncase normalVector::%s():\n{\n', name);
    need = incoming_need_from_static(solveShears);
    print_incoming_second_order_calc(need);
    print_all_moments_block(sol, phase_field);
    fprintf('\n   return;\n}\n');

end

function [Q, w, cx, cy, cz, bitLists, bcs] = defineStencil(stencil)

    if stencil == "D3Q19"
        Q = 19;
        w = zeros(Q, 1);
        w(1) = 1/3; % i=0
        w(2:7) = 1/18; % i=1..6
        w(8:19) = 1/36; % i=7..18

        cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0]';
        cy = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1]';
        cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1]';

        bitLists = { ...
                        [0, 2, 4, 6, 8, 10, 12], ...
                        [0, 1, 4, 6, 12, 13, 15], ...
                        [0, 2, 3, 6, 10, 14, 17], ...
                        [0, 1, 3, 6, 7, 15, 17], ...
                        [0, 2, 4, 5, 8, 16, 18], ...
                        [0, 1, 4, 5, 9, 13, 18], ...
                        [0, 2, 3, 5, 11, 14, 16], ...
                        [0, 1, 3, 5, 7, 9, 11] ...
                    };
    elseif stencil == "D3Q27"
        Q = 27;
        w = zeros(Q, 1);
        w(1) = 8/27; % i=0
        w(2:7) = 2/27; % i=1..6
        w(8:19) = 1/54; % i=7..18
        w(20:27) = 1/216; % i=19..26

        cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1]';
        cy = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1]';
        cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1]';

        bitLists = { ...
                        [0, 2, 4, 6, 8, 10, 12, 20], ...
                        [0, 1, 4, 6, 12, 13, 15, 26], ...
                        [0, 2, 3, 6, 10, 14, 17, 24], ...
                        [0, 1, 3, 6, 7, 15, 17, 21], ...
                        [0, 2, 4, 5, 8, 16, 18, 22], ...
                        [0, 1, 4, 5, 9, 13, 18, 23], ...
                        [0, 2, 3, 5, 11, 14, 16, 25], ...
                        [0, 1, 3, 5, 7, 9, 11, 19] ...
                    };
    else
        error('stencil must be "D3Q19" or "D3Q27".');
    end

    % numeric bitmasks
    bcs = struct( ...
        'SOUTH_WEST_BACK', 128, 'SOUTH_WEST_FRONT', 8, ...
        'SOUTH_EAST_BACK', 64, 'SOUTH_EAST_FRONT', 4, ...
        'SOUTH_WEST', 136, 'SOUTH_EAST', 68, ...
        'WEST_BACK', 160, 'WEST_FRONT', 10, ...
        'EAST_BACK', 80, 'EAST_FRONT', 5, ...
        'SOUTH_BACK', 192, 'SOUTH_FRONT', 12, ...
        'WEST', 170, 'EAST', 85, ...
        'SOUTH', 204, 'BACK', 240, ...
        'FRONT', 15, 'NORTH', 51, ...
        'NORTH_WEST_BACK', 32, 'NORTH_WEST_FRONT', 2, ...
        'NORTH_EAST_BACK', 16, 'NORTH_EAST_FRONT', 1, ...
        'NORTH_BACK', 48, 'NORTH_FRONT', 3, ...
        'NORTH_EAST', 17, 'NORTH_WEST', 34 ...
    );
end

function Os = getCoordinates(n, bitLists)

    if n < 0 || n > 255
        error('Invalid input: n must be between 0 and 255');
    end

    bits = find(bitget(uint16(n), 1:8)) - 1;
    coords = [];

    for k = 1:numel(bits)
        coords = [coords, bitLists{bits(k) + 1}]; %#ok<AGROW>
    end

    Os = unique(sort(coords));
end

function Is = getOppositeCoordinates(Os, cx, cy, cz)
    dirs = [cx(:), cy(:), cz(:)];
    sel = dirs(Os + 1, :);
    opp = -sel;

    Is = [];

    for k = 1:size(opp, 1)
        m = find(dirs(:, 1) == opp(k, 1) & dirs(:, 2) == opp(k, 2) & dirs(:, 3) == opp(k, 3), 1, 'first');

        if ~isempty(m)
            Is(end + 1) = m - 1; %#ok<AGROW>
        end

    end

    Is = unique(sort(Is));
end

function sol = solve_static_case(Q, wS, Hxx, Hxy, Hxz, Hyy, Hyz, Hzz, Os, Is, as2, as4, omega, pI, mxyI, mxzI, myzI, solveShears)
    % Unknowns
    p = sym('p', 'real');

    solveShears = string(solveShears);
    solveShears = solveShears(:).';

    mxy = sym(0); mxz = sym(0); myz = sym(0);
    if any(solveShears == "mxy"); mxy = sym('mxy', 'real'); end
    if any(solveShears == "mxz"); mxz = sym('mxz', 'real'); end
    if any(solveShears == "myz"); myz = sym('myz', 'real'); end

    % static reset
    ux = sym(0); uy = sym(0); uz = sym(0);
    mxx = sym(0); myy = sym(0); mzz = sym(0);

    % build f and feq
    f = sym(zeros(Q, 1));
    feq = sym(zeros(Q, 1));

    for i = 0:(Q - 1)
        ii = i + 1;

        f(ii) = wS(ii) * ( ...
            p + ...
            as2 * (ux * 0 + uy * 0 + uz * 0) + ...
            (as4 / 2) * mxx * Hxx(ii) + as4 * mxy * Hxy(ii) + as4 * mxz * Hxz(ii) + ...
            (as4 / 2) * myy * Hyy(ii) + as4 * myz * Hyz(ii) + (as4 / 2) * mzz * Hzz(ii));

        feq(ii) = wS(ii) * ( ...
            p + ...
            as2 * (ux * 0 + uy * 0 + uz * 0) + ...
            (as4 / 2) * (ux * ux) * Hxx(ii) + as4 * (ux * uy) * Hxy(ii) + as4 * (ux * uz) * Hxz(ii) + ...
            (as4 / 2) * (uy * uy) * Hyy(ii) + as4 * (uy * uz) * Hyz(ii) + (as4 / 2) * (uz * uz) * Hzz(ii));
    end

    % equations
    eqs = sym([]);

    eqPress = pI == sum((1 - omega) * f(Os + 1) + omega * feq(Os + 1));
    eqs(end + 1) = eqPress;

    unknowns = p;

    if any(solveShears == "mxy")
        eqMxy = mxyI == sum(f(Is + 1) .* Hxy(Is + 1));
        eqs(end + 1) = eqMxy;
        unknowns(end + 1) = sym('mxy', 'real');
    end

    if any(solveShears == "mxz")
        eqMxz = mxzI == sum(f(Is + 1) .* Hxz(Is + 1));
        eqs(end + 1) = eqMxz;
        unknowns(end + 1) = sym('mxz', 'real');
    end

    if any(solveShears == "myz")
        eqMyz = myzI == sum(f(Is + 1) .* Hyz(Is + 1));
        eqs(end + 1) = eqMyz;
        unknowns(end + 1) = sym('myz', 'real');
    end

    S = solve(eqs, unknowns);

    sol = struct();

    if isstruct(S)

        if ~isfield(S, 'p') || isempty(S.p)
            error('solve_static_case: no solution for p');
        end

        sol.p = simplify(S.p(1));

        if any(solveShears == "mxy")
            sol.mxy = simplify(S.mxy(1));
        end

        if any(solveShears == "mxz")
            sol.mxz = simplify(S.mxz(1));
        end

        if any(solveShears == "myz")
            sol.myz = simplify(S.myz(1));
        end

    else
        S = sym(S);

        nVars = numel(unknowns);

        if isvector(S)
            Srow = reshape(S, 1, []);
        else

            if size(S, 2) == nVars
                Srow = S(1, :);
            elseif size(S, 1) == nVars
                Srow = transpose(S(:, 1));
            else
                error('solve_static_case: unexpected solve() output size [%d x %d] for nVars=%d', size(S, 1), size(S, 2), nVars);
            end

        end

        sol.p = simplify(Srow(1));

        k = 2;

        if any(solveShears == "mxy")
            sol.mxy = simplify(Srow(k)); k = k + 1;
        end

        if any(solveShears == "mxz")
            sol.mxz = simplify(Srow(k)); k = k + 1;
        end

        if any(solveShears == "myz")
            sol.myz = simplify(Srow(k)); k = k + 1;
        end

    end

end

function s = cpp_sum_pop(indices0)
    parts = arrayfun(@(i) sprintf('pop[q_i<%d>()]', i), indices0, 'UniformOutput', false);
    s = strjoin(parts, ' + ');

    if isempty(s)
        s = 'static_cast<scalar_t>(0)';
    end

end

function s = cpp_sum_pop_weighted(indices0, weights)
    out = "";
    first = true;

    for k = 1:numel(indices0)
        c = weights(k);
        c = double(vpa(c));
        if c == 0, continue; end
        term = sprintf('pop[q_i<%d>()]', indices0(k));

        if first

            if c == 1
                out = string(term);
            elseif c == -1
                out = "-" + string(term);
            else
                out = cpp_cast_num(sym(c)) + " * " + string(term);
            end

            first = false;
        else

            if c == 1
                out = out + " + " + string(term);
            elseif c == -1
                out = out + " - " + string(term);
            else

                if c > 0
                    out = out + " + " + cpp_cast_num(sym(c)) + " * " + string(term);
                else
                    out = out + " - " + cpp_cast_num(sym(-c)) + " * " + string(term);
                end

            end

        end

    end

    if strlength(out) == 0
        s = 'static_cast<scalar_t>(0)';
    else
        s = char(out);
    end

end

function s = cpp_expr(expr)
    expr = simplify(expr);

    if isequal(expr, sym(0))
        s = 'static_cast<scalar_t>(0)';
        return;
    end

    vars = [sym('pI'), sym('mxyI'), sym('mxzI'), sym('myzI'), sym('omega')];

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

    numStr = cpp_poly(Nint, vars);
    denStr = cpp_poly(Dint, vars);

    if isequal(Dint, sym(1))
        s = numStr;
    else
        s = ['(' numStr ') / (' denStr ')'];
    end

end

function s = cpp_poly(P, vars)
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
        mon = cpp_replace_vars(mon);
        mon = strrep(mon, '*', ' * ');

        if strcmp(mon, '1')
            term = string(cpp_cast_num(ak));
        else

            if isequal(ak, sym(1))
                term = string(mon);
            else
                term = string(cpp_cast_num(ak)) + " * " + string(mon);
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

function s = cpp_cast_num(x)
    [xn, xd] = numden(x);

    if isequal(xd, sym(1))
        s = sprintf('static_cast<scalar_t>(%s)', char(xn));
    else
        s = sprintf('(static_cast<scalar_t>(%s) / static_cast<scalar_t>(%s))', char(xn), char(xd));
    end

end

function s = cpp_replace_vars(s)
    if isstring(s), s = char(s); end

    s = regexprep(s, '\<ux\>\^2', 'ux*ux');
    s = regexprep(s, '\<uy\>\^2', 'uy*uy');
    s = regexprep(s, '\<uz\>\^2', 'uz*uz');

    s = regexprep(s, '\<myzI\>', 'myz_I');
    s = regexprep(s, '\<mxzI\>', 'mxz_I');
    s = regexprep(s, '\<mxyI\>', 'mxy_I');
    s = regexprep(s, '\<mzzI\>', 'mzz_I');
    s = regexprep(s, '\<myyI\>', 'myy_I');
    s = regexprep(s, '\<mxxI\>', 'mxx_I');
    s = regexprep(s, '\<pI\>', 'p_I');

    s = regexprep(s, '\<omega\>', 'device::omega');

    s = regexprep(s, '\<p\>', 'moments[m_i<0>()]');
    s = regexprep(s, '\<ux\>', 'moments[m_i<1>()]');
    s = regexprep(s, '\<uy\>', 'moments[m_i<2>()]');
    s = regexprep(s, '\<uz\>', 'moments[m_i<3>()]');
end

function print_all_moments_block(sol, phase_field)

    mxy = sym(0); if isfield(sol, 'mxy'); mxy = sol.mxy; end
    mxz = sym(0); if isfield(sol, 'mxz'); mxz = sol.mxz; end
    myz = sym(0); if isfield(sol, 'myz'); myz = sol.myz; end

    fprintf('   moments[m_i<0>()] = %s; // p\n', cpp_expr(sol.p));
    fprintf('   moments[m_i<1>()] = static_cast<scalar_t>(0); // ux\n');
    fprintf('   moments[m_i<2>()] = static_cast<scalar_t>(0); // uy\n');
    fprintf('   moments[m_i<3>()] = static_cast<scalar_t>(0); // uz\n');
    fprintf('   moments[m_i<4>()] = static_cast<scalar_t>(0); // mxx\n');
    fprintf('   moments[m_i<5>()] = %s; // mxy\n', cpp_expr(mxy));
    fprintf('   moments[m_i<6>()] = %s; // mxz\n', cpp_expr(mxz));
    fprintf('   moments[m_i<7>()] = static_cast<scalar_t>(0); // myy\n');
    fprintf('   moments[m_i<8>()] = %s; // myz\n', cpp_expr(myz));
    fprintf('   moments[m_i<9>()] = static_cast<scalar_t>(0); // mzz\n');

    if phase_field == true
        fprintf('   moments[m_i<10>()] = static_cast<scalar_t>(0); // phi\n');
    end

end

function need = incoming_need_init()
    need = struct('mxx', false, 'myy', false, 'mzz', false, ...
        'mxy', false, 'mxz', false, 'myz', false);
end

function print_incoming_second_order_calc(need)
    % Prints only the required 2nd-order incoming moments.
    % Intentionally does NOT print p_I.

    if ~(need.mxx || need.myy || need.mzz || need.mxy || need.mxz || need.myz)
        return;
    end

    fprintf('   // Incoming moments\n');

    if need.mxx
        fprintf('   const scalar_t mxx_I = velocitySet::calculate_moment<VelocitySet, X, X>(pop, boundaryNormal);\n');
    end

    if need.myy
        fprintf('   const scalar_t myy_I = velocitySet::calculate_moment<VelocitySet, Y, Y>(pop, boundaryNormal);\n');
    end

    if need.mzz
        fprintf('   const scalar_t mzz_I = velocitySet::calculate_moment<VelocitySet, Z, Z>(pop, boundaryNormal);\n');
    end

    if need.mxy
        fprintf('   const scalar_t mxy_I = velocitySet::calculate_moment<VelocitySet, X, Y>(pop, boundaryNormal);\n');
    end

    if need.mxz
        fprintf('   const scalar_t mxz_I = velocitySet::calculate_moment<VelocitySet, X, Z>(pop, boundaryNormal);\n');
    end

    if need.myz
        fprintf('   const scalar_t myz_I = velocitySet::calculate_moment<VelocitySet, Y, Z>(pop, boundaryNormal);\n');
    end

    fprintf('\n');
end

function need = incoming_need_from_static(solveShears)
    need = incoming_need_init();

    solveShears = string(solveShears);
    solveShears = solveShears(:).';

    % Static generator only ever uses shear incoming moments.
    need.mxy = any(solveShears == "mxy");
    need.mxz = any(solveShears == "mxz");
    need.myz = any(solveShears == "myz");
end
