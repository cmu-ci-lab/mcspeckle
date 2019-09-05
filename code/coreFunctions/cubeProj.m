function [xd]=cubeProj(x,vmin,vmax,v)

    thr=0.1^7;
    bd = inf(size(v));
    for j=1:size(vmin,1)
        bd(j,:)=(x(j)-vmax(j))./max(abs(v(j,:)),thr).*(v(j,:)<=0)+(-x(j)+vmin(j))./max(abs(v(j,:)),thr).*(v(j,:)>0);
    end

    d = max(bd);

    if d>=0
        ep=0.1^4*max(abs(vmax-vmin));

        xd=x+(d+ep)*v;
    else
        xd=d;
    end

end