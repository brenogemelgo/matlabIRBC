function [dx, dy, dz] = neumannShift(name)
    dx = 0; dy = 0; dz = 0;

    if contains(name, "NORTH")
        dy = -1;
    end

    if contains(name, "SOUTH")
        dy = +1;
    end

    if contains(name, "WEST")
        dx = +1;
    end

    if contains(name, "EAST")
        dx = -1;
    end

    if contains(name, "BACK")
        dz = +1;
    end

    if contains(name, "FRONT")
        dz = -1;
    end

end
