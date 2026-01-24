function incomingSecondOrder(need)

    if ~(need.mxx || need.myy || need.mzz || need.mxy || need.mxz || need.myz)
        return;
    end

    fprintf('   // Incoming moments\n');

    if need.mxx
        fprintf('   const scalar_t mxx_I = velocitySet::calculate_moment<VelocitySet, axis::X, axis::X>(pop, boundaryNormal) * inv_rho_I;\n');
    end

    if need.myy
        fprintf('   const scalar_t myy_I = velocitySet::calculate_moment<VelocitySet, axis::Y, axis::Y>(pop, boundaryNormal) * inv_rho_I;\n');
    end

    if need.mzz
        fprintf('   const scalar_t mzz_I = velocitySet::calculate_moment<VelocitySet, axis::Z, axis::Z>(pop, boundaryNormal) * inv_rho_I;\n');
    end

    if need.mxy
        fprintf('   const scalar_t mxy_I = velocitySet::calculate_moment<VelocitySet, axis::X, axis::Y>(pop, boundaryNormal) * inv_rho_I;\n');
    end

    if need.mxz
        fprintf('   const scalar_t mxz_I = velocitySet::calculate_moment<VelocitySet, axis::X, axis::Z>(pop, boundaryNormal) * inv_rho_I;\n');
    end

    if need.myz
        fprintf('   const scalar_t myz_I = velocitySet::calculate_moment<VelocitySet, axis::Y, axis::Z>(pop, boundaryNormal) * inv_rho_I;\n');
    end

    fprintf('\n');
end
