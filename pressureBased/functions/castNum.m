function s = castNum(x)
    [xn, xd] = numden(x);

    if isequal(xd, sym(1))
        s = sprintf('static_cast<scalar_t>(%s)', char(xn));
    else
        s = sprintf('(static_cast<scalar_t>(%s) / static_cast<scalar_t>(%s))', char(xn), char(xd));
    end

end
