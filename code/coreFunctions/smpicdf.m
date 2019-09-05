function x=smpicdf(icdf)

    N=length(icdf);

    a=rand;
    i=ceil(a*N);
    x=icdf(i);

end