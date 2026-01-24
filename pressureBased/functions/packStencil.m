function S = packStencil(stencil)
    [Q, w, cx, cy, cz, bitLists, bcs] = defineStencil(stencil);

    S = struct();
    S.Q = Q;
    S.cx = cx(:); S.cy = cy(:); S.cz = cz(:);
    S.cxS = sym(S.cx); S.cyS = sym(S.cy); S.czS = sym(S.cz);
    S.wS = sym(w(:));
    S.bitLists = bitLists;
    S.bcs = bcs;

    S.as = sqrt(sym(3));
    S.as2 = S.as ^ 2;
    S.as4 = S.as ^ 4;
    S.cs2 = sym(1) / S.as2;

    S.Hxx = S.cxS .^ 2 - S.cs2;
    S.Hxy = S.cxS .* S.cyS;
    S.Hxz = S.cxS .* S.czS;
    S.Hyy = S.cyS .^ 2 - S.cs2;
    S.Hyz = S.cyS .* S.czS;
    S.Hzz = S.czS .^ 2 - S.cs2;
end
