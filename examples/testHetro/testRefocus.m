
lambda=1;
box_min=[-50; -50];
box_max=[ 50;  50];
sigt=[1/25,1/25,0,1/50;0,1/100,1/100,1/25]';
box_bin=[25;50];

theta_l=pi-deg2rad([-35:5:35]);
theta_v=deg2rad([-40:0.2:40]);


Nl=length(theta_l);
Nv=length(theta_v);


Nx=(box_max(1)-box_min(1))/box_bin(1);
Nz=(box_max(2)-box_min(2))/box_bin(2);
Nxz=Nx*Nz;

v=[sin(theta_v);cos(theta_v)];
l=[sin(theta_l);cos(theta_l)];

doCBS=0; smpFlg=2;
maxItr=10^4;
%maxItr=1;
g=0.5
j=1
rng(520)
tic
uL(:,:,j)=MCsampleHetro( sigt, 1, box_min,box_max,box_bin, l, v,1,1,maxItr,lambda,doCBS,smpFlg,1,0);
toc
box_stp=[0.5;1];
for j=1:Nl
    tic
    rU(:,:,j)=refocus(uL(:,j),v,l(:,j),box_min,box_max,box_stp,lambda);
    toc
    figure, imshow(abs(rU(:,:,j)),[])
    %keyboard
end

figure, imshow(abs(mean(rU,3)),[])