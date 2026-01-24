function printMoment(idx, tag, expr19, expr27)

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
