clc; clearvars; close all

root = fileparts(mfilename("fullpath"));
addpath(fullfile(root, "functions"));

phase_field = true;

wantNames = [ ...
                 "WEST_SOUTH_BACK", "WEST_SOUTH_FRONT", "EAST_SOUTH_BACK", "EAST_SOUTH_FRONT", ...
                 "WEST_NORTH_BACK", "WEST_NORTH_FRONT", "EAST_NORTH_BACK", "EAST_NORTH_FRONT", ...
                 ...
                 "WEST_SOUTH", "EAST_SOUTH", "WEST_NORTH", "EAST_NORTH", ...
                 "WEST_BACK", "WEST_FRONT", "EAST_BACK", "EAST_FRONT", ...
                 "SOUTH_BACK", "SOUTH_FRONT", "NORTH_BACK", "NORTH_FRONT", ...
                 ...
                 "WEST", "EAST", "SOUTH", "NORTH", "BACK", "FRONT"
             ];

S19 = buildStencilData("D3Q19");
S27 = buildStencilData("D3Q27");

% symbols
pI = sym('pI', 'real');
mxyI = sym('mxyI', 'real');
mxzI = sym('mxzI', 'real');
myzI = sym('myzI', 'real');
omega = sym('omega', 'real');

% cases
corner = { ...
              struct('name', "WEST_SOUTH_BACK", 'key', "WEST_SOUTH_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "WEST_SOUTH_FRONT", 'key', "WEST_SOUTH_FRONT", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_SOUTH_BACK", 'key', "EAST_SOUTH_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_SOUTH_FRONT", 'key', "EAST_SOUTH_FRONT", 'solve', string.empty(1, 0)), ...
              struct('name', "WEST_NORTH_BACK", 'key', "WEST_NORTH_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "WEST_NORTH_FRONT", 'key', "WEST_NORTH_FRONT", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_NORTH_BACK", 'key', "EAST_NORTH_BACK", 'solve', string.empty(1, 0)), ...
              struct('name', "EAST_NORTH_FRONT", 'key', "EAST_NORTH_FRONT", 'solve', string.empty(1, 0)) ...
          };

edge = { ...
            struct('name', "WEST_SOUTH", 'key', "WEST_SOUTH", 'solve', ["mxy"]), ...
            struct('name', "EAST_SOUTH", 'key', "EAST_SOUTH", 'solve', ["mxy"]), ...
            struct('name', "WEST_NORTH", 'key', "WEST_NORTH", 'solve', ["mxy"]), ...
            struct('name', "EAST_NORTH", 'key', "EAST_NORTH", 'solve', ["mxy"]), ...
            ...
            struct('name', "WEST_BACK", 'key', "WEST_BACK", 'solve', ["mxz"]), ...
            struct('name', "WEST_FRONT", 'key', "WEST_FRONT", 'solve', ["mxz"]), ...
            struct('name', "EAST_BACK", 'key', "EAST_BACK", 'solve', ["mxz"]), ...
            struct('name', "EAST_FRONT", 'key', "EAST_FRONT", 'solve', ["mxz"]), ...
            ...
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

mask = cellfun(@(s) any(s.name == wantNames), allCases);
allCases = allCases(mask);

% main loop
for c = 1:numel(allCases)
    name = allCases{c}.name;
    key = allCases{c}.key;

    solveShears = string(allCases{c}.solve);
    solveShears = solveShears(:).';

    Os19 = getCoordinates(S19.bcs.(key), S19.bitLists);
    Is19 = getOppositeCoordinates(Os19, S19.cx, S19.cy, S19.cz);

    Os27 = getCoordinates(S27.bcs.(key), S27.bitLists);
    Is27 = getOppositeCoordinates(Os27, S27.cx, S27.cy, S27.cz);

    sol19 = solve_static_case(S19.Q, S19.wS, S19.Hxx, S19.Hxy, S19.Hxz, S19.Hyy, S19.Hyz, S19.Hzz, ...
        Os19, Is19, S19.as2, S19.as4, omega, pI, mxyI, mxzI, myzI, solveShears);

    sol27 = solve_static_case(S27.Q, S27.wS, S27.Hxx, S27.Hxy, S27.Hxz, S27.Hyy, S27.Hyz, S27.Hzz, ...
        Os27, Is27, S27.as2, S27.as4, omega, pI, mxyI, mxzI, myzI, solveShears);

    fprintf('\ncase normalVector::%s():\n{\n', name);

    need = incoming_need_from_static(solveShears);
    print_incoming_second_order_calc(need);

    print_all_moments_block_conditional(sol19, sol27, phase_field);

    fprintf('\n   return;\n}\n');
end

function sol = solve_static_case(Q, wS, Hxx, Hxy, Hxz, Hyy, Hyz, Hzz, Os, Is, as2, as4, omega, pI, mxyI, mxzI, myzI, solveShears)
    % unknowns
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

function print_all_moments_block_conditional(sol19, sol27, phase_field)

    mxy19 = sym(0); if isfield(sol19, 'mxy'); mxy19 = sol19.mxy; end
    mxz19 = sym(0); if isfield(sol19, 'mxz'); mxz19 = sol19.mxz; end
    myz19 = sym(0); if isfield(sol19, 'myz'); myz19 = sol19.myz; end

    mxy27 = sym(0); if isfield(sol27, 'mxy'); mxy27 = sol27.mxy; end
    mxz27 = sym(0); if isfield(sol27, 'mxz'); mxz27 = sol27.mxz; end
    myz27 = sym(0); if isfield(sol27, 'myz'); myz27 = sol27.myz; end

    print_moment_conditional(0, 'p', sol19.p, sol27.p);

    print_moment_conditional(1, 'ux', sym(0), sym(0));
    print_moment_conditional(2, 'uy', sym(0), sym(0));
    print_moment_conditional(3, 'uz', sym(0), sym(0));

    print_moment_conditional(4, 'mxx', sym(0), sym(0));

    print_moment_conditional(5, 'mxy', mxy19, mxy27);
    print_moment_conditional(6, 'mxz', mxz19, mxz27);

    print_moment_conditional(7, 'myy', sym(0), sym(0));
    print_moment_conditional(8, 'myz', myz19, myz27);
    print_moment_conditional(9, 'mzz', sym(0), sym(0));

    if phase_field == true
        print_moment_conditional(10, 'phi', sym(0), sym(0));
    end

end

function print_moment_conditional(idx, tag, expr19, expr27)

    if exprEqual(expr19, expr27)
        fprintf('   moments[m_i<%d>()] = %s; // %s\n', idx, expression(expr19), tag);
        return;
    end

    fprintf('   if constexpr (VelocitySet::Q() == 19)\n');
    fprintf('   {\n');
    fprintf('       moments[m_i<%d>()] = %s; // %s\n', idx, expression(expr19), tag);
    fprintf('   }\n');
    fprintf('   else\n');
    fprintf('   {\n');
    fprintf('       moments[m_i<%d>()] = %s; // %s\n', idx, expression(expr27), tag);
    fprintf('   }\n');
end

function need = incoming_need_init()
    need = struct('mxx', false, 'myy', false, 'mzz', false, ...
        'mxy', false, 'mxz', false, 'myz', false);
end

function print_incoming_second_order_calc(need)

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

    need.mxy = any(solveShears == "mxy");
    need.mxz = any(solveShears == "mxz");
    need.myz = any(solveShears == "myz");
end
