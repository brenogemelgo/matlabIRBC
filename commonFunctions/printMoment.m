function printMoment(idx, tag, expr19, expr27, mode)

    if nargin < 5
        mode = "auto";
    end

    mode = string(mode);

    if mode == "19"
        fprintf('       moments[m_i<%d>()] = %s; // %s\n', idx, expression(expr19), tag);
        return;
    elseif mode == "27"
        fprintf('       moments[m_i<%d>()] = %s; // %s\n', idx, expression(expr27), tag);
        return;
    elseif mode ~= "auto"
        error('printMoment: unknown mode "%s". Use "auto", "19", or "27".', mode);
    end

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
