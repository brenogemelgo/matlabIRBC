function e = getOr(sol, field, fallback)
    if isfield(sol, field), e = sol.(field); else, e = fallback; end
end
