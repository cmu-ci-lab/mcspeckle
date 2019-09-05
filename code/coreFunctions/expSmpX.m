function [x,px]=expSmpX(box_min,box_max,lmean,sigs)
    % exponential distribution position of first scatter
    dim=size(lmean,1); x=zeros(dim,1);   
    lz=abs(lmean(end));
    sigs=sigs/lz;
    dmax=(box_max(end)-box_min(end));
    rmax=1-exp(-dmax*sigs);
    d=-log(-rand.*rmax+1)/(sigs);
    x(1:end-1)=rand(size(box_min(1:end-1))).*(box_max(1:end-1)-box_min(1:end-1))+box_min(1:end-1);
    x(end)=(d+box_min(end)).*(lmean(end)>0)+(-d+box_max(end)).*(lmean(end)<=0);

    px = exp(-sigs*d)/(rmax)*sigs*dmax;
end