
a=deg2rad(-90);
v=-[sin(a); cos(a)];
x=[-50;-30];
lambda=1
box_min=[-60;-60]; box_max=[60;60]; box_bin=[40;60];
sigt=[1/200,0,1/400;0,1/400,1/200]';

is_ff_v=1;
ff_sign=-1


s=[0:1:500];
e_v=zeros(1,length(s));
e_v0=(evalphaseattHetro(x,ff_sign*v,is_ff_v,sigt,lambda,box_min,box_max,box_bin))^2;
ts=sigt(findBin(x,box_min,box_max,box_bin))
ep=0.1^6;
for j=1:length(s)
    if(max((x+s(j)*v)>(box_max-ep))|max(x+s(j)*v<(box_min+ep)))
          break
    end
    ts=sigt(findBin(x+s(j)*v,box_min,box_max,box_bin));
    e_v(j)=ts*e_v0./(evalphaseattHetro(x+s(j)*v,ff_sign*v,is_ff_v,sigt,lambda,box_min,box_max,box_bin))^2;
end

for j=1:10^6
    d2(j)=smpWoodcock(x,box_min,box_max,box_bin,sigt,v);
end

h2=hist(d2,s+0.5);
h2=h2/sum(h2);

figure, plot(s, abs(e_v))
hold on
plot(s,h2/(h2(1))*max(e_v(1)))
