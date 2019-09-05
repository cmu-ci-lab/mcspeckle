function funcVals = evaluateHG(theta, gVec, cosineFlag, dim)
    % HG scattering rotation amplitude
    if(nargin < 2)
        gVec = 0;
    end

    if(size(gVec, 2) ~= 1)
        gVec = gVec';
    end

    if(nargin < 3)
        cosineFlag = 0;
    end

    if(nargin < 4)
        dim = 3;
    end

    if(cosineFlag == 0)
        cosTheta = cos(theta);
    else
        cosTheta = theta;
    end

    if(dim == 3)
        funcVals = 1 / 4 / pi *...
            bsxfun(@rdivide,...
            1 - gVec .^ 2, (1 + bsxfun(@plus,...
            gVec .^ 2, - 2 * bsxfun(@times,...
            gVec, cosTheta))) .^ (3 / 2));
    end

    if(dim == 2)
        funcVals = 1 / 2 / pi *...
            bsxfun(@rdivide,...
            1 - gVec .^ 2, (1 + bsxfun(@plus,...
            gVec .^ 2, - 2 * bsxfun(@times,...
            gVec, cosTheta))));
    end

end