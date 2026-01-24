function Os = getCoordinates(n, bitLists)

    if n < 0 || n > 255
        error('Invalid input: n must be between 0 and 255');
    end

    bits = find(bitget(uint16(n), 1:8)) - 1;
    coords = [];

    for k = 1:numel(bits)
        coords = [coords, bitLists{bits(k) + 1}]; %#ok<AGROW>
    end

    Os = unique(sort(coords));
end
