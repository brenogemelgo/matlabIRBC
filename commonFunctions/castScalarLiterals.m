function s = castScalarLiterals(s)
    % Wrap numeric literals with static_cast<scalar_t>(...)
    % - avoids double-wrapping numbers already inside static_cast<scalar_t>(...)
    % - avoids touching tokens adjacent to < > (template args) or identifiers

    pattern = '(?<!scalar_t>\()(?<![A-Za-z0-9_<>])(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(?![A-Za-z0-9_<>])';
    s = regexprep(s, pattern, 'static_cast<scalar_t>($1)');
end
