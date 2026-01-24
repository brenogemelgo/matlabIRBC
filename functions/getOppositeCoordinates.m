function Is = getOppositeCoordinates(Os, cx, cy, cz)
    dirs = [cx(:), cy(:), cz(:)];
    sel = dirs(Os + 1, :);
    opp = -sel;

    Is = [];

    for k = 1:size(opp, 1)
        m = find(dirs(:, 1) == opp(k, 1) & dirs(:, 2) == opp(k, 2) & dirs(:, 3) == opp(k, 3), 1, 'first');

        if ~isempty(m)
            Is(end + 1) = m - 1; %#ok<AGROW>
        end

    end

    Is = unique(sort(Is));
end
