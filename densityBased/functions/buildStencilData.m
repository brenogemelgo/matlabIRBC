function S = buildStencilData(stencil)
    [Q, w, cx, cy, cz, bitLists, bcs] = defineStencil(stencil);

    S = struct();
    S.Q = Q;
    S.cx = cx(:); S.cy = cy(:); S.cz = cz(:);
    S.bitLists = bitLists;
    S.bcs = bcs;

    S.as = sqrt(sym(3));
    S.as2 = S.as ^ 2;
    S.as4 = S.as ^ 4;

    cxS = sym(S.cx); cyS = sym(S.cy); czS = sym(S.cz);
    S.wS = sym(w(:));

    S.Hxx = cxS .^ 2 - S.as ^ (-2);
    S.Hxy = cxS .* cyS;
    S.Hxz = cxS .* czS;
    S.Hyy = cyS .^ 2 - S.as ^ (-2);
    S.Hyz = cyS .* czS;
    S.Hzz = czS .^ 2 - S.as ^ (-2);
end
