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

S19 = packStencil("D3Q19");
S27 = packStencil("D3Q27");

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
              struct('name', "WEST_SOUTH_BACK", 'key', "WEST_SOUTH_BACK", 'type', "corner"), ...
              struct('name', "WEST_SOUTH_FRONT", 'key', "WEST_SOUTH_FRONT", 'type', "corner"), ...
              struct('name', "EAST_SOUTH_BACK", 'key', "EAST_SOUTH_BACK", 'type', "corner"), ...
              struct('name', "EAST_SOUTH_FRONT", 'key', "EAST_SOUTH_FRONT", 'type', "corner"), ...
              struct('name', "WEST_NORTH_BACK", 'key', "WEST_NORTH_BACK", 'type', "corner"), ...
              struct('name', "WEST_NORTH_FRONT", 'key', "WEST_NORTH_FRONT", 'type', "corner"), ...
              struct('name', "EAST_NORTH_BACK", 'key', "EAST_NORTH_BACK", 'type', "corner"), ...
              struct('name', "EAST_NORTH_FRONT", 'key', "EAST_NORTH_FRONT", 'type', "corner") ...
          };

edge_xz_yz = { ...
                  struct('name', "WEST_SOUTH", 'key', "WEST_SOUTH", 'type', "edge_xz_yz"), ...
                  struct('name', "EAST_SOUTH", 'key', "EAST_SOUTH", 'type', "edge_xz_yz"), ...
                  struct('name', "WEST_NORTH", 'key', "WEST_NORTH", 'type', "edge_xz_yz"), ...
                  struct('name', "EAST_NORTH", 'key', "EAST_NORTH", 'type', "edge_xz_yz") ...
              };

edge_xy_yz = { ...
                  struct('name', "WEST_BACK", 'key', "WEST_BACK", 'type', "edge_xy_yz"), ...
                  struct('name', "WEST_FRONT", 'key', "WEST_FRONT", 'type', "edge_xy_yz"), ...
                  struct('name', "EAST_BACK", 'key', "EAST_BACK", 'type', "edge_xy_yz"), ...
                  struct('name', "EAST_FRONT", 'key', "EAST_FRONT", 'type', "edge_xy_yz") ...
              };

edge_xy_xz = { ...
                  struct('name', "SOUTH_BACK", 'key', "SOUTH_BACK", 'type', "edge_xy_xz"), ...
                  struct('name', "SOUTH_FRONT", 'key', "SOUTH_FRONT", 'type', "edge_xy_xz"), ...
                  struct('name', "NORTH_BACK", 'key', "NORTH_BACK", 'type', "edge_xy_xz"), ...
                  struct('name', "NORTH_FRONT", 'key', "NORTH_FRONT", 'type', "edge_xy_xz") ...
              };

face_we = { ...
               struct('name', "WEST", 'key', "WEST", 'type', "face_we"), ...
               struct('name', "EAST", 'key', "EAST", 'type', "face_we") ...
           };

face_sn = { ...
               struct('name', "SOUTH", 'key', "SOUTH", 'type', "face_sn"), ...
               struct('name', "NORTH", 'key', "NORTH", 'type', "face_sn") ...
           };

face_bf = { ...
               struct('name', "BACK", 'key', "BACK", 'type', "face_bf"), ...
               struct('name', "FRONT", 'key', "FRONT", 'type', "face_bf") ...
           };

allCases = [corner(:); edge_xz_yz(:); edge_xy_yz(:); edge_xy_xz(:); face_we(:); face_sn(:); face_bf(:)];
allCases = allCases.';

mask = cellfun(@(s) any(s.name == wantNames), allCases);
allCases = allCases(mask);

% main loop
for c = 1:numel(allCases)
    name = allCases{c}.name;
    key = allCases{c}.key;
    typ = allCases{c}.type;

    Os19 = getCoordinates(S19.bcs.(key), S19.bitLists);
    Is19 = getOppositeCoordinates(Os19, S19.cx, S19.cy, S19.cz);

    Os27 = getCoordinates(S27.bcs.(key), S27.bitLists);
    Is27 = getOppositeCoordinates(Os27, S27.cx, S27.cy, S27.cz);

    sol19 = solve_specific_case(S19.Q, S19.wS, S19.Hxx, S19.Hxy, S19.Hxz, S19.Hyy, S19.Hyz, S19.Hzz, ...
        Os19, Is19, S19.as2, S19.as4, S19.cxS, S19.cyS, S19.czS, ...
        p, ux, uy, uz, pI, mxxI, myyI, mzzI, mxyI, mxzI, myzI, typ);

    sol27 = solve_specific_case(S27.Q, S27.wS, S27.Hxx, S27.Hxy, S27.Hxz, S27.Hyy, S27.Hyz, S27.Hzz, ...
        Os27, Is27, S27.as2, S27.as4, S27.cxS, S27.cyS, S27.czS, ...
        p, ux, uy, uz, pI, mxxI, myyI, mzzI, mxyI, mxzI, myzI, typ);

    print_neumann_case_block_conditional(name, typ, sol19, sol27, phase_field);

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
        case "edge_xz_yz"
            mxz = sym('mxz', 'real'); myz = sym('myz', 'real');
            unknowns = [mxz; myz];

        case "edge_xy_yz"
            mxy = sym('mxy', 'real'); myz = sym('myz', 'real');
            unknowns = [mxy; myz];

        case "edge_xy_xz"
            mxy = sym('mxy', 'real'); mxz = sym('mxz', 'real');
            unknowns = [mxy; mxz];

        case "face_we"
            myy = sym('myy', 'real'); mxy = sym('mxy', 'real'); mxz = sym('mxz', 'real');
            myz = sym('myz', 'real'); mzz = sym('mzz', 'real');
            unknowns = [myy; mxy; mxz; myz; mzz];

        case "face_sn"
            mxx = sym('mxx', 'real'); mxy = sym('mxy', 'real'); mxz = sym('mxz', 'real');
            myz = sym('myz', 'real'); mzz = sym('mzz', 'real');
            unknowns = [mxx; mxy; mxz; myz; mzz];

        case "face_bf"
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
        case "edge_xz_yz"
            eqs = [ ...
                       mxzI == sumHxz; ...
                       myzI == sumHyz ...
                   ];

        case "edge_xy_yz"
            eqs = [ ...
                       mxyI == sumHxy; ...
                       myzI == sumHyz ...
                   ];

        case "edge_xy_xz"
            eqs = [ ...
                       mxyI == sumHxy; ...
                       mxzI == sumHxz ...
                   ];

        case "face_we"
            eqmYY = (myy == uy * uy + uz * uz - mzz);
            eqmZZ = (myyI - mzzI == sumHyy - sumHzz);
            eqXY = (mxyI == sumHxy);
            eqXZ = (mxzI == sumHxz);
            eqYZ = (myzI == sumHyz);

            eqs = [eqmYY; eqXY; eqXZ; eqYZ; eqmZZ];

        case "face_sn"
            eqmXX = (mxx == ux * ux + uz * uz - mzz);
            eqmZZ = (mxxI - mzzI == sumHxx - sumHzz);
            eqXY = (mxyI == sumHxy);
            eqXZ = (mxzI == sumHxz);
            eqYZ = (myzI == sumHyz);

            eqs = [eqmXX; eqXY; eqXZ; eqYZ; eqmZZ];

        case "face_bf"
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

function print_neumann_case_block_conditional(name, typ, sol19, sol27, phase_field)

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

    fprintf('\n   // IRBC-Neumann\n');

    ux = sym('ux', 'real'); uy = sym('uy', 'real'); uz = sym('uz', 'real');

    e19 = moment_exprs(sol19, ux, uy, uz);
    e27 = moment_exprs(sol27, ux, uy, uz);

    print_moment_conditional(4, 'mxx', e19.mxx, e27.mxx);
    print_moment_conditional(5, 'mxy', e19.mxy, e27.mxy);
    print_moment_conditional(6, 'mxz', e19.mxz, e27.mxz);
    print_moment_conditional(7, 'myy', e19.myy, e27.myy);
    print_moment_conditional(8, 'myz', e19.myz, e27.myz);
    print_moment_conditional(9, 'mzz', e19.mzz, e27.mzz);

    fprintf('\n   return;\n');
    fprintf('}\n');
end

function E = moment_exprs(sol, ux, uy, uz)
    E.mxx = get_or(sol, 'mxx', ux * ux);
    E.mxy = get_or(sol, 'mxy', ux * uy);
    E.mxz = get_or(sol, 'mxz', ux * uz);
    E.myy = get_or(sol, 'myy', uy * uy);
    E.myz = get_or(sol, 'myz', uy * uz);
    E.mzz = get_or(sol, 'mzz', uz * uz);
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
        case "edge_xz_yz"
            need.mxz = true;
            need.myz = true;

        case "edge_xy_yz"
            need.mxy = true;
            need.myz = true;

        case "edge_xy_xz"
            need.mxy = true;
            need.mxz = true;

        case "face_we"
            % uses myy_I and mzz_I plus all 3 shears
            need.myy = true;
            need.mzz = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "face_sn"
            % uses mxx_I and mzz_I plus all 3 shears
            need.mxx = true;
            need.mzz = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "face_bf"
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
