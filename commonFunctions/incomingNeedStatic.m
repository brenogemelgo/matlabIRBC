function need = incomingNeedStatic(solveShears)
    need = incomingNeedInit();

    solveShears = string(solveShears);
    solveShears = solveShears(:).';

    need.mxy = any(solveShears == "mxy");
    need.mxz = any(solveShears == "mxz");
    need.myz = any(solveShears == "myz");
end
