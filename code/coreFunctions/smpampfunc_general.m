function w = smpampfunc_general(ow, sct_type,ampfunc)
    % sample new scatteing direction

    dim = size(ow,1);
    switch sct_type
        case 1
        w = randn(dim,1);
        w = w/norm(w);
        
        case 2
        a = smpampfunc(ampfunc);
        cosa = cos(a);

        case 3
        cosa = sampleHG(ampfunc,dim);
    end

    if sct_type>1
        w = rotateBytheta(ow,cosa);
    end

end