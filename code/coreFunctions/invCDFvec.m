function icdf=invCDFvec(cdf,N)

    if ~exist('N','var')
        N=length(cdf);
    end

    icdf =cumsum(hist(cdf,N));
end