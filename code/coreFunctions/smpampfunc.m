function ang=smpampfunc(ampfunc)
    % sample direction by provided amplitude function
    N=length(ampfunc.sampleIcdf);
    stp=1/(N-1);
    a=rand;
    i=round(a/stp)+1;
    ang=ampfunc.sampleIcdf(i);
end
