function need = incomingNeedNeumann(typ)
    need = incomingNeedInit();

    switch typ
        case "edge_xz_yz"
            need.mxz = true;
            need.myz = true;

        case "edge_xy_yz"
            need.mxy = true;
            need.myz = true;

        case "edge_xy_xz"
            need.mxy = true;
            need.mxz = true;

        case "face_we"
            need.myy = true;
            need.mzz = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "face_sn"
            need.mxx = true;
            need.mzz = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "face_bf"
            need.mxx = true;
            need.myy = true;
            need.mxy = true;
            need.mxz = true;
            need.myz = true;

        case "corner"
            % nothing

        otherwise
            error('incomingNeedNeumann: unknown typ "%s".', typ);
    end

end
