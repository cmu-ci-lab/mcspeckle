maxItr=10^6;

sigt=[1/30,0,1/60;0,1/120,1/30]';
box_bin=[40;60]/5;
sigt=imresize(sigt,5,'nearest');

H=zeros(size(sigt));
P=zeros(size(sigt));
lmean=[0;1];
lmean=mean(l,2);
lmean=[0.3;0.8]
lmean=lmean/norm(lmean)
lrad=sqrt(sum((box_max-box_min).^2));

for itr=1:maxItr
   [x,px]=expSmpXHetro(box_min,box_max,box_bin,lmean,sigt,lrad);
                        
   
   if(max(x>box_max)|max(x<box_min))
        continue
   end
 
   
    bi=findBin(x,box_min,box_max,box_bin);
 
    H(bi)=H(bi)+1;
end
    
H=H/maxItr/prod(box_bin);


stp=0.5;
[gx,gz]=ndgrid([box_min(1)+stp/2:stp:box_max(1)],[box_min(2)+stp/2:stp:box_max(2)]);

for itr=1:length(gx(:))
    x=[gx(itr);gz(itr)];
    px=abs(evalphaseattHetro(x,lmean,1,sigt,1,box_min,box_max,box_bin))^2;

    px=px/(lrad)^(1);

    bi=findBin(x,box_min,box_max,box_bin);
    P(bi)=P(bi)+px;
end
P=P/prod(box_bin/stp);
P=P.*sigt;
P
H

sum(P)*10
 sum(H)*10