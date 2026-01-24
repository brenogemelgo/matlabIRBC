clc; clearvars; close all

stencil = "D3Q27"; % "D3Q19" or "D3Q27"
phase_field = true;

wantNames = [ ...
                 "WEST_SOUTH_BACK", "WEST_SOUTH_FRONT", "EAST_SOUTH_BACK", "EAST_SOUTH_FRONT", "WEST_NORTH_BACK", "WEST_NORTH_FRONT", "EAST_NORTH_BACK", "EAST_NORTH_FRONT", ...
                 "SOUTH_WEST", "SOUTH_EAST", "NORTH_WEST", "NORTH_EAST", "WEST_BACK", "WEST_FRONT", "EAST_BACK", "EAST_FRONT", "SOUTH_BACK", "SOUTH_FRONT", "NORTH_BACK", "NORTH_FRONT", ...
                 "WEST", "EAST", "SOUTH", "NORTH", "BACK", "FRONT" ...
             ];

[Q, w, cx, cy, cz, bitLists, bcs] = defineStencil(stencil);

% constants
as = sqrt(sym(3));
as2 = as ^ 2;
as4 = as ^ 4;
cs2 = sym(1) / as2; % = 1/3
cxS = sym(cx(:));
cyS = sym(cy(:));
czS = sym(cz(:));
wS = sym(w(:));

% hermite tensors
Hxx = cxS .^ 2 - cs2;
Hxy = cxS .* cyS;
Hxz = cxS .* czS;
Hyy = cyS .^ 2 - cs2;
Hyz = cyS .* czS;
Hzz = czS .^ 2 - cs2;

% symbols
p = sym('p', 'real');
ux = sym('ux', 'real');
uy = sym('uy', 'real');
uz = sym('uz', 'real');
pI = sym('pI', 'real');
mxxI = sym('mxxI', 'real');
myyI = sym('myyI', 'real');
mzzI = sym('mzzI', 'real');
mxyI = sym('mxyI', 'real');
mxzI = sym('mxzI', 'real');
myzI = sym('myzI', 'real');

% cases
corner = { ...
              struct('name', "WEST_SOUTH_BACK", 'key', "SOUTH_WEST_BACK", 'type', "corner"), ...
              struct('name', "WEST_SOUTH_FRONT", 'key', "SOUTH_WEST_FRONT", 'type', "corner"), ...
              struct('name', "EAST_SOUTH_BACK", 'key', "SOUTH_EAST_BACK", 'type', "corner"), ...
              struct('name', "EAST_SOUTH_FRONT", 'key', "SOUTH_EAST_FRONT", 'type', "corner"), ...
              struct('name', "WEST_NORTH_BACK", 'key', "NORTH_WEST_BACK", 'type', "corner"), ...
              struct('name', "WEST_NORTH_FRONT", 'key', "NORTH_WEST_FRONT", 'type', "corner"), ...
              struct('name', "EAST_NORTH_BACK", 'key', "NORTH_EAST_BACK", 'type', "corner"), ...
              struct('name', "EAST_NORTH_FRONT", 'key', "NORTH_EAST_FRONT", 'type', "corner") ...
          };

edge_mxz_myz = { ...
                    struct('name', "SOUTH_WEST", 'key', "SOUTH_WEST", 'type', "edge_mxz_myz"), ...
                    struct('name', "SOUTH_EAST", 'key', "SOUTH_EAST", 'type', "edge_mxz_myz"), ...
                    struct('name', "NORTH_WEST", 'key', "NORTH_WEST", 'type', "edge_mxz_myz"), ...
                    struct('name', "NORTH_EAST", 'key', "NORTH_EAST", 'type', "edge_mxz_myz") ...
                };

edge_mxy_myz = { ...
                    struct('name', "WEST_BACK", 'key', "WEST_BACK", 'type', "edge_mxy_myz"), ...
                    struct('name', "WEST_FRONT", 'key', "WEST_FRONT", 'type', "edge_mxy_myz"), ...
                    struct('name', "EAST_BACK", 'key', "EAST_BACK", 'type', "edge_mxy_myz"), ...
                    struct('name', "EAST_FRONT", 'key', "EAST_FRONT", 'type', "edge_mxy_myz") ...
                };

edge_mxy_mxz = { ...
                    struct('name', "SOUTH_BACK", 'key', "SOUTH_BACK", 'type', "edge_mxy_mxz"), ...
                    struct('name', "SOUTH_FRONT", 'key', "SOUTH_FRONT", 'type', "edge_mxy_mxz"), ...
                    struct('name', "NORTH_BACK", 'key', "NORTH_BACK", 'type', "edge_mxy_mxz"), ...
                    struct('name', "NORTH_FRONT", 'key', "NORTH_FRONT", 'type', "edge_mxy_mxz") ...
                };

face_west_east = { ...
                      struct('name', "WEST", 'key', "WEST", 'type', "face_west_east"), ...
                      struct('name', "EAST", 'key', "EAST", 'type', "face_west_east") ...
                  };

face_south_north = { ...
                        struct('name', "SOUTH", 'key', "SOUTH", 'type', "face_south_north"), ...
                        struct('name', "NORTH", 'key', "NORTH", 'type', "face_south_north") ...
                    };

face_back_front = { ...
                       struct('name', "BACK", 'key', "BACK", 'type', "face_back_front"), ...
                       struct('name', "FRONT", 'key', "FRONT", 'type', "face_back_front") ...
                   };

allCases = [corner(:); edge_mxz_myz(:); edge_mxy_myz(:); edge_mxy_mxz(:); face_west_east(:); face_south_north(:); face_back_front(:)];
allCases = allCases.';

mask = cellfun(@(s) any(s.name == wantNames), allCases);
allCases = allCases(mask);

% main loop
for c = 1:numel(allCases)
    name = allCases{c}.name;
    key = allCases{c}.key;
    typ = allCases{c}.type;

    Os = getCoordinates(bcs.(key), bitLists);
    Is = getOppositeCoordinates(Os, cx, cy, cz);

    % solve boundary moments
    sol = solve_specific_case(Q, wS, Hxx, Hxy, Hxz, Hyy, Hyz, Hzz, Os, Is, ...
        as2, as4, cxS, cyS, czS, p, ux, uy, uz, pI, mxxI, myyI, mzzI, mxyI, mxzI, myzI, typ);

    % print solved boundary moments
    print_neumann_case_block(name, typ, Is, cs2, cx, cy, cz, Hxy, Hxz, Hyz, sol, phase_field);

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

function sol = solve_specific_case(Q, wS, Hxx, Hxy, Hxz, Hyy, Hyz, Hzz, Os, Is, as2, as4, ...
        cxS, cyS, czS, ...
        p, ux, uy, uz, ...
        pI, mxxI, myyI, mzzI, mxyI, mxzI, myzI, typ) %#ok<INUSD>

    % equilibrium moments for anything not explicitly solved
    mxx = ux * ux;
    myy = uy * uy;
    mzz = uz * uz;
    mxy = ux * uy;
    mxz = ux * uz;
    myz = uy * uz;

    % corner: equilibrium
    if typ == "corner"
        sol = struct();
        return;
    end

    % promote only the specified moments to symbolic unknowns
    switch typ
        case "edge_mxz_myz"
            mxz = sym('mxz', 'real'); myz = sym('myz', 'real');
            unknowns = [mxz; myz];

        case "edge_mxy_myz"
            mxy = sym('mxy', 'real'); myz = sym('myz', 'real');
            unknowns = [mxy; myz];

        case "edge_mxy_mxz"
            mxy = sym('mxy', 'real'); mxz = sym('mxz', 'real');
            unknowns = [mxy; mxz];

        case "face_west_east"
            myy = sym('myy', 'real'); mxy = sym('mxy', 'real'); mxz = sym('mxz', 'real');
            myz = sym('myz', 'real'); mzz = sym('mzz', 'real');
            unknowns = [myy; mxy; mxz; myz; mzz];

        case "face_south_north"
            mxx = sym('mxx', 'real'); mxy = sym('mxy', 'real'); mxz = sym('mxz', 'real');
            myz = sym('myz', 'real'); mzz = sym('mzz', 'real');
            unknowns = [mxx; mxy; mxz; myz; mzz];

        case "face_back_front"
            mxx = sym('mxx', 'real'); mxy = sym('mxy', 'real'); mxz = sym('mxz', 'real');
            myy = sym('myy', 'real'); myz = sym('myz', 'real');
            unknowns = [mxx; mxy; mxz; myy; myz];

        otherwise
            error('solve_specific_case: unknown typ "%s".', typ);
    end

    % build f
    f = sym(zeros(Q, 1));

    for i = 0:(Q - 1)
        ii = i + 1;

        f(ii) = wS(ii) * ( ...
            p + ...
            as2 * (ux * cxS(ii) + uy * cyS(ii) + uz * czS(ii)) + ...
            (as4 / 2) * mxx * Hxx(ii) + as4 * mxy * Hxy(ii) + as4 * mxz * Hxz(ii) + ...
            (as4 / 2) * myy * Hyy(ii) + as4 * myz * Hyz(ii) + (as4 / 2) * mzz * Hzz(ii));
    end

    idx = Is + 1; %
    sumHxy = sum(f(idx) .* Hxy(idx));
    sumHxz = sum(f(idx) .* Hxz(idx));
    sumHyz = sum(f(idx) .* Hyz(idx));
    sumHxx = sum(f(idx) .* Hxx(idx));
    sumHyy = sum(f(idx) .* Hyy(idx));
    sumHzz = sum(f(idx) .* Hzz(idx));

    % build equations and unknowns
    switch typ
        case "edge_mxz_myz"
            eqs = [ ...
                       mxzI == sumHxz; ...
                       myzI == sumHyz ...
                   ];

        case "edge_mxy_myz"
            eqs = [ ...
                       mxyI == sumHxy; ...
                       myzI == sumHyz ...
                   ];

        case "edge_mxy_mxz"
            eqs = [ ...
                       mxyI == sumHxy; ...
                       mxzI == sumHxz ...
                   ];

        case "face_west_east"
            eqmYY = (myy == uy * uy + uz * uz - mzz);
            eqmZZ = (myyI - mzzI == sumHyy - sumHzz);
            eqXY = (mxyI == sumHxy);
            eqXZ = (mxzI == sumHxz);
            eqYZ = (myzI == sumHyz);

            eqs = [eqmYY; eqXY; eqXZ; eqYZ; eqmZZ];

        case "face_south_north"
            eqmXX = (mxx == ux * ux + uz * uz - mzz);
            eqmZZ = (mxxI - mzzI == sumHxx - sumHzz);
            eqXY = (mxyI == sumHxy);
            eqXZ = (mxzI == sumHxz);
            eqYZ = (myzI == sumHyz);

            eqs = [eqmXX; eqXY; eqXZ; eqYZ; eqmZZ];

        case "face_back_front"
            eqmXX = (mxx == ux * ux + uy * uy - myy);
            eqmYY = (mxxI - myyI == sumHxx - sumHyy);
            eqXY = (mxyI == sumHxy);
            eqXZ = (mxzI == sumHxz);
            eqYZ = (myzI == sumHyz);

            eqs = [eqmXX; eqXY; eqXZ; eqmYY; eqYZ];

        otherwise
            error('solve_specific_case: unknown typ "%s".', typ);
    end

    S = solve(eqs, unknowns);

    sol = struct();

    if isstruct(S)

        for k = 1:numel(unknowns)
            vn = char(unknowns(k));

            if ~isfield(S, vn) || isempty(S.(vn))
                error('solve_specific_case: no solution for "%s" (typ=%s).', vn, typ);
            end

            sol.(vn) = simplify(S.(vn)(1));
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
                error('solve_specific_case: unexpected solve() output size [%d x %d] for nVars=%d (typ=%s).', ...
                    size(S, 1), size(S, 2), nVars, typ);
            end

        end

        for k = 1:nVars
            vn = char(unknowns(k));
            sol.(vn) = simplify(Srow(k));
        end

    end

end

function s = cpp_sum_pop(indices0)
    parts = arrayfun(@(i) sprintf('pop[q_i<%d>()]', i), indices0, 'UniformOutput', false);
    s = strjoin(parts, ' + ');
    if isempty(s), s = 'static_cast<scalar_t>(0)'; end
end

function s = cpp_sum_pop_weighted(indices0, weights)
    out = "";
    first = true;

    for k = 1:numel(indices0)
        c = weights(k);
        if isa(c, 'sym'), c = double(vpa(c)); else, c = double(c); end

        if abs(c - round(c)) < 1e-12
            c = round(c);
        end

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

    vars = [ ...
                sym('p'), sym('ux'), sym('uy'), sym('uz'), ...
                sym('pI'), sym('mxxI'), sym('myyI'), sym('mzzI'), ...
                sym('mxyI'), sym('mxzI'), sym('myzI') ...
            ];

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
            if neg, out = "-" + term; else, out = term; end
            first = false;
        else
            if neg, out = out + " - " + term; else, out = out + " + " + term; end
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

function print_neumann_case_block(name, typ, Is, cs2, cx, cy, cz, Hxy, Hxz, Hyz, sol, phase_field)

    [dx, dy, dz] = neumann_shift_from_name(name);

    fprintf('\ncase normalVector::%s():\n{\n', name);
    fprintf('   const label_t tid = device::idxBlock(threadIdx.x%s, threadIdx.y%s, threadIdx.z%s);\n\n', ...
        off_str(dx), off_str(dy), off_str(dz));

    fprintf('   // Classic Neumann\n');
    fprintf('   moments[m_i<0>()] = shared_buffer[tid * (NUMBER_MOMENTS<true>() + 1) + m_i<0>()];   // p\n');
    fprintf('   moments[m_i<1>()] = shared_buffer[tid * (NUMBER_MOMENTS<true>() + 1) + m_i<1>()];   // ux\n');
    fprintf('   moments[m_i<2>()] = shared_buffer[tid * (NUMBER_MOMENTS<true>() + 1) + m_i<2>()];   // uy\n');
    fprintf('   moments[m_i<3>()] = shared_buffer[tid * (NUMBER_MOMENTS<true>() + 1) + m_i<3>()];   // uz\n');

    if phase_field
        fprintf('   moments[m_i<10>()] = shared_buffer[tid * (NUMBER_MOMENTS<true>() + 1) + m_i<10>()]; // phi\n');
    end

    need = incoming_need_from_neumann_type(typ);
    print_incoming_second_order_calc(need);

    fprintf('\n');

    % explicit incoming moments
    % print_incoming_block_neumann(typ, Is, cs2, cx, cy, cz, Hxy, Hxz, Hyz);

    fprintf('   // IRBC-Neumann\n');

    % unspecified moments default to equilibrium u_alpha u_beta
    ux = sym('ux', 'real'); uy = sym('uy', 'real'); uz = sym('uz', 'real');

    expr_mxx = get_or(sol, 'mxx', ux * ux);
    expr_mxy = get_or(sol, 'mxy', ux * uy);
    expr_mxz = get_or(sol, 'mxz', ux * uz);
    expr_myy = get_or(sol, 'myy', uy * uy);
    expr_myz = get_or(sol, 'myz', uy * uz);
    expr_mzz = get_or(sol, 'mzz', uz * uz);

    fprintf('   moments[m_i<4>()] = %s; // mxx\n', cpp_expr(expr_mxx));
    fprintf('   moments[m_i<5>()] = %s; // mxy\n', cpp_expr(expr_mxy));
    fprintf('   moments[m_i<6>()] = %s; // mxz\n', cpp_expr(expr_mxz));
    fprintf('   moments[m_i<7>()] = %s; // myy\n', cpp_expr(expr_myy));
    fprintf('   moments[m_i<8>()] = %s; // myz\n', cpp_expr(expr_myz));
    fprintf('   moments[m_i<9>()] = %s; // mzz\n\n', cpp_expr(expr_mzz));

    fprintf('   return;\n');
    fprintf('}\n');
end

function print_incoming_block_neumann(typ, Is, cs2, cx, cy, cz, Hxy, Hxz, Hyz)
    fprintf('   // Incoming moments\n');

    switch typ
        case "edge_mxz_myz"
            fprintf('   const scalar_t mxz_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxz(Is + 1))));
            fprintf('   const scalar_t myz_I = %s;\n\n', cpp_sum_pop_weighted(Is, double(Hyz(Is + 1))));

        case "edge_mxy_myz"
            fprintf('   const scalar_t mxy_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxy(Is + 1))));
            fprintf('   const scalar_t myz_I = %s;\n\n', cpp_sum_pop_weighted(Is, double(Hyz(Is + 1))));

        case "edge_mxy_mxz"
            fprintf('   const scalar_t mxy_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxy(Is + 1))));
            fprintf('   const scalar_t mxz_I = %s;\n\n', cpp_sum_pop_weighted(Is, double(Hxz(Is + 1))));

        case "face_west_east"
            fprintf('   const scalar_t mxy_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxy(Is + 1))));
            fprintf('   const scalar_t mxz_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxz(Is + 1))));
            fprintf('   const scalar_t myz_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hyz(Is + 1))));
            fprintf('   const scalar_t myy_I = %s - %s * moments[m_i<0>()];\n', ...
                cpp_sum_pop_weighted(Is, double(cy(Is + 1) .^ 2)), cpp_cast_num(cs2));
            fprintf('   const scalar_t mzz_I = %s - %s * moments[m_i<0>()];\n\n', ...
                cpp_sum_pop_weighted(Is, double(cz(Is + 1) .^ 2)), cpp_cast_num(cs2));

        case "face_south_north"
            fprintf('   const scalar_t mxy_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxy(Is + 1))));
            fprintf('   const scalar_t mxz_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxz(Is + 1))));
            fprintf('   const scalar_t myz_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hyz(Is + 1))));
            fprintf('   const scalar_t mxx_I = %s - %s * moments[m_i<0>()];\n', ...
                cpp_sum_pop_weighted(Is, double(cx(Is + 1) .^ 2)), cpp_cast_num(cs2));
            fprintf('   const scalar_t mzz_I = %s - %s * moments[m_i<0>()];\n\n', ...
                cpp_sum_pop_weighted(Is, double(cz(Is + 1) .^ 2)), cpp_cast_num(cs2));

        case "face_back_front"
            fprintf('   const scalar_t mxy_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxy(Is + 1))));
            fprintf('   const scalar_t mxz_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hxz(Is + 1))));
            fprintf('   const scalar_t myz_I = %s;\n', cpp_sum_pop_weighted(Is, double(Hyz(Is + 1))));
            fprintf('   const scalar_t mxx_I = %s - %s * moments[m_i<0>()];\n', ...
                cpp_sum_pop_weighted(Is, double(cx(Is + 1) .^ 2)), cpp_cast_num(cs2));
            fprintf('   const scalar_t myy_I = %s - %s * moments[m_i<0>()];\n\n', ...
                cpp_sum_pop_weighted(Is, double(cy(Is + 1) .^ 2)), cpp_cast_num(cs2));

        otherwise
            fprintf('\n');
    end

end

function e = get_or(sol, field, fallback)
    if isfield(sol, field), e = sol.(field); else, e = fallback; end
end

function s = off_str(d)

    if d == 0
        s = '';
    elseif d > 0
        s = sprintf(' + %d', d);
    else
        s = sprintf(' - %d', -d);
    end

end

function [dx, dy, dz] = neumann_shift_from_name(name)
    dx = 0; dy = 0; dz = 0;

    if contains(name, "NORTH")
        dy = -1;
    end

    if contains(name, "SOUTH")
        dy = +1;
    end

    if contains(name, "WEST")
        dx = +1;
    end

    if contains(name, "EAST")
        dx = -1;
    end

    if contains(name, "BACK")
        dz = +1;
    end

    if contains(name, "FRONT")
        dz = -1;
    end

end

function need = incoming_need_init()
    need = struct('mxx', false, 'myy', false, 'mzz', false, ...
        'mxy', false, 'mxz', false, 'myz', false);
end

function print_incoming_second_order_calc(need)

    if ~(need.mxx || need.myy || need.mzz || need.mxy || need.mxz || need.myz)
        return;
    end

    fprintf('\n');
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

end

function need = incoming_need_from_neumann_type(typ)
    need = incoming_need_init();

    switch typ
        case "edge_mxz_myz"
            need.mxz = true;
            need.myz = true;

        case "edge_mxy_myz"
            need.mxy = true;
            need.myz = true;

        case "edge_mxy_mxz"
            need.mxy = true;
            need.mxz = true;

        case "face_west_east"
            % uses myy_I and mzz_I plus all 3 shears
            need.myy = true;
            need.mzz = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "face_south_north"
            % uses mxx_I and mzz_I plus all 3 shears
            need.mxx = true;
            need.mzz = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "face_back_front"
            % uses mxx_I and myy_I plus all 3 shears
            need.mxx = true;
            need.myy = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "corner"
            % nothing

        otherwise
            error('incoming_need_from_neumann_type: unknown typ "%s".', typ);
    end

end
