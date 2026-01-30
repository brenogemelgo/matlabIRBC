clc; clearvars; close all

root = fileparts(mfilename("fullpath"));
addpath(fullfile(root, "functions"));
addpath(fullfile(root, "..", "commonFunctions"));

phase_field = true;

% every
% wantNames = [ ...
%                  "WEST_SOUTH_BACK", "WEST_SOUTH_FRONT", "EAST_SOUTH_BACK", "EAST_SOUTH_FRONT", ...
%                  "WEST_NORTH_BACK", "WEST_NORTH_FRONT", "EAST_NORTH_BACK", "EAST_NORTH_FRONT", ...
%                  ...
%                  "WEST_SOUTH", "EAST_SOUTH", "WEST_NORTH", "EAST_NORTH", ...
%                  "WEST_BACK", "WEST_FRONT", "EAST_BACK", "EAST_FRONT", ...
%                  "SOUTH_BACK", "SOUTH_FRONT", "NORTH_BACK", "NORTH_FRONT", ...
%                  ...
%                  "WEST", "EAST", "SOUTH", "NORTH", "BACK", "FRONT"
%              ];

% jet
% wantNames = [ ...
%                  "WEST_SOUTH_FRONT", "WEST_NORTH_FRONT", ...
%                  "EAST_SOUTH_FRONT", "EAST_NORTH_FRONT", ...
%                  ...
%                  "WEST_FRONT", "EAST_FRONT", ...
%                  "SOUTH_FRONT", "NORTH_FRONT", ...
%                  ...
%                  "FRONT"
%              ];

% ssmd
wantNames = [ ...
                 "WEST_NORTH_FRONT", "EAST_NORTH_FRONT", ...
                 ...
                 "WEST_NORTH", "EAST_NORTH", ...
                 "WEST_FRONT", "EAST_FRONT", ...
                 "NORTH_FRONT", ...
                 ...
                 "NORTH", "FRONT"
             ];

S19 = packStencil("D3Q19");
S27 = packStencil("D3Q27");

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

for c = 1:numel(allCases)
    name = allCases{c}.name;
    key = allCases{c}.key;
    typ = allCases{c}.type;

    Os19 = getCoordinates(S19.bcs.(key), S19.bitLists);
    Is19 = getOppositeCoordinates(Os19, S19.cx, S19.cy, S19.cz);

    Os27 = getCoordinates(S27.bcs.(key), S27.bitLists);
    Is27 = getOppositeCoordinates(Os27, S27.cx, S27.cy, S27.cz);

    sol19 = solveNeumann(S19.Q, S19.wS, S19.Hxx, S19.Hxy, S19.Hxz, S19.Hyy, S19.Hyz, S19.Hzz, ...
        Os19, Is19, S19.as2, S19.as4, S19.cxS, S19.cyS, S19.czS, ...
        p, ux, uy, uz, pI, mxxI, myyI, mzzI, mxyI, mxzI, myzI, typ);

    sol27 = solveNeumann(S27.Q, S27.wS, S27.Hxx, S27.Hxy, S27.Hxz, S27.Hyy, S27.Hyz, S27.Hzz, ...
        Os27, Is27, S27.as2, S27.as4, S27.cxS, S27.cyS, S27.czS, ...
        p, ux, uy, uz, pI, mxxI, myyI, mzzI, mxyI, mxzI, myzI, typ);

    printAllMoments(name, typ, sol19, sol27, phase_field);

end

function sol = solveNeumann(Q, wS, Hxx, Hxy, Hxz, Hyy, Hyz, Hzz, Os, Is, as2, as4, ...
        cxS, cyS, czS, ...
        p, ux, uy, uz, ...
        pI, mxxI, myyI, mzzI, mxyI, mxzI, myzI, typ) %#ok<INUSD>

    mxx = ux * ux;
    myy = uy * uy;
    mzz = uz * uz;
    mxy = ux * uy;
    mxz = ux * uz;
    myz = uy * uz;

    if typ == "corner"
        sol = struct();
        return;
    end

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
            error('solveNeumann: unknown typ "%s".', typ);
    end

    f = sym(zeros(Q, 1));

    for i = 0:(Q - 1)
        ii = i + 1;

        f(ii) = wS(ii) * (p + ...
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
            error('solveNeumann: unknown typ "%s".', typ);
    end

    S = solve(eqs, unknowns);

    sol = struct();

    if isstruct(S)

        for k = 1:numel(unknowns)
            vn = char(unknowns(k));

            if ~isfield(S, vn) || isempty(S.(vn))
                error('solveNeumann: no solution for "%s" (typ=%s).', vn, typ);
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
                error('solveNeumann: unexpected solve() output size [%d x %d] for nVars=%d (typ=%s).', ...
                    size(S, 1), size(S, 2), nVars, typ);
            end

        end

        for k = 1:nVars
            vn = char(unknowns(k));
            sol.(vn) = simplify(Srow(k));
        end

    end

end

function printAllMoments(name, typ, sol19, sol27, phase_field)

    [dx, dy, dz] = neumannShift(name);

    fprintf('\ncase normalVector::%s():\n{\n', name);
    fprintf('   const label_t tid = device::idxBlock(threadIdx.x%s, threadIdx.y%s, threadIdx.z%s);\n\n', ...
        offStr(dx), offStr(dy), offStr(dz));

    if phase_field
        nm = 'NUMBER_MOMENTS<true>()';
    else
        nm = 'NUMBER_MOMENTS<false>()';
    end

    fprintf('   // Classic Neumann\n');
    fprintf('   moments[m_i<0>()] = shared_buffer[tid * (%s + 1) + m_i<0>()];   // p\n', nm);
    fprintf('   moments[m_i<1>()] = shared_buffer[tid * (%s + 1) + m_i<1>()];   // ux\n', nm);
    fprintf('   moments[m_i<2>()] = shared_buffer[tid * (%s + 1) + m_i<2>()];   // uy\n', nm);
    fprintf('   moments[m_i<3>()] = shared_buffer[tid * (%s + 1) + m_i<3>()];   // uz\n', nm);

    if phase_field
        fprintf('   moments[m_i<10>()] = shared_buffer[tid * (NUMBER_MOMENTS<true>() + 1) + m_i<10>()]; // phi\n');
    end

    need = incomingNeedNeumann(typ);
    incoming2ndNeumann(need);

    fprintf('\n   // IRBC-Neumann\n');

    ux = sym('ux', 'real'); uy = sym('uy', 'real'); uz = sym('uz', 'real');

    e19 = momentExpr(sol19, ux, uy, uz);
    e27 = momentExpr(sol27, ux, uy, uz);

    same2 = ...
        exprEqual(e19.mxx, e27.mxx) && ...
        exprEqual(e19.mxy, e27.mxy) && ...
        exprEqual(e19.mxz, e27.mxz) && ...
        exprEqual(e19.myy, e27.myy) && ...
        exprEqual(e19.myz, e27.myz) && ...
        exprEqual(e19.mzz, e27.mzz);

    if same2
        printMoment(4, 'mxx', e19.mxx, e27.mxx);
        printMoment(5, 'mxy', e19.mxy, e27.mxy);
        printMoment(6, 'mxz', e19.mxz, e27.mxz);
        printMoment(7, 'myy', e19.myy, e27.myy);
        printMoment(8, 'myz', e19.myz, e27.myz);
        printMoment(9, 'mzz', e19.mzz, e27.mzz);
    else
        fprintf('   if constexpr (VelocitySet::Q() == 19)\n');
        fprintf('   {\n');
        printMoment(4, 'mxx', e19.mxx, e27.mxx, "19");
        printMoment(5, 'mxy', e19.mxy, e27.mxy, "19");
        printMoment(6, 'mxz', e19.mxz, e27.mxz, "19");
        printMoment(7, 'myy', e19.myy, e27.myy, "19");
        printMoment(8, 'myz', e19.myz, e27.myz, "19");
        printMoment(9, 'mzz', e19.mzz, e27.mzz, "19");
        fprintf('   }\n');
        fprintf('   else\n');
        fprintf('   {\n');
        printMoment(4, 'mxx', e19.mxx, e27.mxx, "27");
        printMoment(5, 'mxy', e19.mxy, e27.mxy, "27");
        printMoment(6, 'mxz', e19.mxz, e27.mxz, "27");
        printMoment(7, 'myy', e19.myy, e27.myy, "27");
        printMoment(8, 'myz', e19.myz, e27.myz, "27");
        printMoment(9, 'mzz', e19.mzz, e27.mzz, "27");
        fprintf('   }\n');
    end

    fprintf('\n   return;\n');
    fprintf('}\n');
end
