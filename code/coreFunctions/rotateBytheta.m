function w=rotateBytheta(ow,costheta)
    % get new vector direction by given theta direction. phi direction is
    % sampled uniformly
    sintheta = real(sqrt(1-real(costheta)^2));
    dim = size(ow,1);
    if(dim == 2)
        sintheta = ((rand>0.5)*2-1)*sintheta;
        w = [costheta,sintheta;-sintheta,costheta]*ow;
    end

    if(dim == 3)
        phi = pi*2*rand;

        sinphi = sin(phi);
        cosphi = cos(phi);
        w = zeros(3,1);
        if (abs(ow(3)) > 0.999999999)
            w(1) = sintheta * cosphi;
            w(2) = sintheta * sinphi;
            w(3) = costheta * ow(3) / abs(ow(3));
        else
            temp = sqrt(1 - ow(3)^2);
            w(1) = sintheta*(ow(1)*ow(3)*cosphi - ow(2)*sinphi)/temp + ow(1)*costheta;
            w(2) = sintheta*(ow(2)*ow(3)*cosphi + ow(1)*sinphi)/temp + ow(2)*costheta;
            w(3) = -sintheta*cosphi*temp+ow(3)*costheta;
        end

    end
end