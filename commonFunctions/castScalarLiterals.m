function s = castScalarLiterals(s)
    pattern = '(?<!scalar_t>\()(?<![A-Za-z0-9_<>])(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(?![A-Za-z0-9_<>])';
    s = regexprep(s, pattern, 'static_cast<scalar_t>($1)');
end
