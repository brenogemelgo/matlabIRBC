function s = replaceVars(s)
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
    s = regexprep(s, '\<rhoI\>', 'rho_I');

    s = regexprep(s, '\<omega\>', 'device::omega');

    s = regexprep(s, '\<rho\>', 'moments[m_i<0>()]');
    s = regexprep(s, '\<ux\>', 'moments[m_i<1>()]');
    s = regexprep(s, '\<uy\>', 'moments[m_i<2>()]');
    s = regexprep(s, '\<uz\>', 'moments[m_i<3>()]');
end
