function e = evalphaseatt(x,v,is_ff,sigt,lambda,box_min,box_max)
    % complex transmission and attenuation term in multiple scattering
    
    Nv = size(v,2);
    dim = size(v,1);
    if (is_ff)
        pv = x'*v;
        rv = ones(1,Nv);
        bdv = cubeDist(x,box_min,box_max,-v);
    else
        d = repmat(x,1,Nv)-v;
        nd = sum(d.^2,1).^0.5;
        v = d./repmat(nd,dim,1);
        pv = nd; 
        rv = 1./(nd+0.01*lambda).^((dim-1)/2);
        bdv = min(cubeDist(x,box_min,box_max,-v),nd);  
    end

     e = exp(1i*pi*2/lambda*pv-sigt/2*bdv).*rv;
 
end