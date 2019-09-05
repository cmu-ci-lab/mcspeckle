function icdf=invCDF(cdf,dim)
    % inverse of cumulative distribution function
    N = length(cdf);
    stp = 1/(N-1);
    grid = 0:stp:1;
    icdf=zeros(size(cdf));
    for j=1:N
        [~,icdf(j)] = min(abs(cdf - grid(j)));
    end
    
    if(dim == 2)
        icdf = icdf/N*2*pi;
    end
    
    if(dim == 3)
        icdf = icdf/N*pi;
    end

end