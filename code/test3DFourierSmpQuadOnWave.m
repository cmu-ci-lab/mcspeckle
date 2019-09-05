
lambda=1;
box_min=[-50; -50;  0];
box_max=[ 50;  50;  100];
box_min=[-4000; -4000;  0];
box_max=[ 4000;  4000;  100];
box_min=[-1000;-1000;0]; box_max=[1000;1000;30];%3000



sigt=1/400;
sigt=1/40000;
sigt=1/2000;
%sigt=1/50;
%theta_l=deg2rad([-35:5:35]);
theta_l=deg2rad([-20:20:20]);

theta_max=deg2rad(40);
theta_stp=deg2rad(0.2);
theta_max=0.5; 
theta_max=0.02; 
theta_stp=theta_max/200;
theta_l=[-theta_max:theta_stp*50:theta_max];


%IMPORTANT: this is a uniform grid of sin(theta) not of theta
%This is important to have a uniform grid for the Fourier transform
theta_v=[-theta_max:theta_stp:theta_max];

Nl0=length(theta_l);
Nv0=length(theta_v);
Nv=Nv0^2;
Nl=Nl0^2;

%Nx=(box_max(1)-box_min(1))/box_bin(1);
%Nz=(box_max(2)-box_min(2))/box_bin(2);
%Nxz=Nx*Nz;

[vx,vy]=ndgrid(theta_v,theta_v);
v=-[vx(:)';vy(:)'; sqrt(1-vx(:).^2-vy(:).^2)'];  %Note again: we assume the grid is a uniform grid of sin(theta) so no sin applyied  
%v=[sin(vx(:))';sin(vy(:))';ones(1,Nv)];
%v=v./repmat(sum(v.^2,1).^0.5,3,1);
vsign=-1;
[lx,ly]=ndgrid(theta_l,theta_l);
l=-[lx(:)';ly(:)'; sqrt(1-lx(:).^2-ly(:).^2)'];  
%[lx,ly]=ndgrid(theta_l,theta_l);
%l=[sin(lx(:))';sin(ly(:))';ones(1,Nl)];
%l=l./repmat(sum(l.^2,1).^0.5,3,1);




doCBS=1; smpFlg=1;
maxItr=10^4*4;
maxItr=10^3
maxItr=3;
%maxItr=100;
maxItr=50;
j=1

rng(520)
for j=1:1
tic
uLf(:,:,:,j)=MCfieldFourierQuad( sigt, 1, box_min,box_max, l, theta_max,theta_stp,maxItr,lambda,doCBS,smpFlg,[],vsign);



toc
end
uLf=sum(uLf,4);
%keyboard


lW=rand(Nl,1).*exp(2*pi*i*rand(Nl,1));
%lW(:)=0; lW(1)=1;
uLfwp=0;
for j=1:Nl
    uLfwp=uLfwp+uLf(:,:,j)*lW(j);
end
rng(520)
tic
[uLfw]=MCfieldFourierQuadOnWave( sigt, 1, box_min,box_max, l, theta_max,theta_stp,maxItr,lambda,doCBS,smpFlg,[0;0;1],lW,vsign);
toc

max(max(abs(uLfw-uLfwp)))

[meanU]=evalMeanUnifCtr(box_min,box_max,l,v,sigt,lambda);
meanlW=reshape(meanU*lW(:),Nv0,Nv0);

figure, imshow(real(uLfw),[]) 
figure, imshow(real(uLfw+meanlW),[]) 



return 
