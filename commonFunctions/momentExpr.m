function E = momentExpr(sol, ux, uy, uz)
    E.mxx = getOr(sol, 'mxx', ux * ux);
    E.mxy = getOr(sol, 'mxy', ux * uy);
    E.mxz = getOr(sol, 'mxz', ux * uz);
    E.myy = getOr(sol, 'myy', uy * uy);
    E.myz = getOr(sol, 'myz', uy * uz);
    E.mzz = getOr(sol, 'mzz', uz * uz);
end
