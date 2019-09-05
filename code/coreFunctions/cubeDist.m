function [d,xd] = cubeDist(x,vmin,vmax,v)
    % distance from v to x inside the box

    thr=0.1^7;
    bd = inf(size(v));
    for j=1:size(vmin,1)
        bd(j,:) = (x(j)-vmin(j))./max(abs(v(j,:)),thr).*(v(j,:)<=0) + ...
            (-x(j)+vmax(j))./max(abs(v(j,:)),thr).*(v(j,:)>0);
    end

    d=min(bd);
    
    if nargout>1
        ep = 0.1^4*max(abs(vmax-vmin));
        xd = x+(d+ep)*v;
    end
    
end