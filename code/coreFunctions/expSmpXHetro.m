function [x,px]=expSmpXHetro(box_min,box_max,box_bin,lmean,sigt,lrad)
    dim=size(lmean,1); x=zeros(dim,1);   

    if isempty(lrad)
        lrad=sqrt(sum((box_max-box_min).^2));
    end

    u=null(lmean');

    x=u*((rand(dim-1,1)-0.5)*lrad)-lmean*max(abs([box_max;box_min]))*dim;

    x=cubeProj(x,box_min,box_max,lmean);

    if(max(x>box_max)||max(x<box_min))
        px=0;
       return     
    end


    d=smpWoodcock(x,box_min,box_max,box_bin,sigt,lmean);

    x=x+d*lmean;
    if(max(x>box_max)||max(x<box_min))
        px=0;
       return     
    end

    px=abs(evalphaseattHetro(x,lmean,1,sigt,1,box_min,box_max,box_bin))^2;

    px=px/(lrad)^(dim-1);

end