lambda=1
box_min=[-60;-60]; box_max=[60;60]; box_bin=[40;60];
%sigt=[1/200,0,1/400;0,1/400,1/200]';
sigt=[1/30,0,1/60;0,1/120,1/30]';

MFP=1./sigt;





theta_l=pi-deg2rad([0,4,12,32]);
theta_v=deg2rad([-88:4:270]);




Nl=length(theta_l);
Nv=length(theta_v);


Nx=(box_max(1)-box_min(1))/box_bin(1);
Nz=(box_max(2)-box_min(2))/box_bin(2);
Nxz=Nx*Nz;


load tmpdata/ampfunc05lmbd1
%eval(sprintf('load ../../resdirHetro/res_itr36_MFP200.mat'))
eval(sprintf('load ../../resdirHetro/res_itr100_MFP30.mat'))

E2=covU;
for j1=1:Nl
    for j2=1:Nl
        covU(:,:,j1,j2)=covU(:,:,j1,j2)-meanU(:,j1)*meanU(:,j2)';
    end
end


v=[sin(theta_v);cos(theta_v)];
l=[sin(theta_l);cos(theta_l)];

doCBS=1; smpFlg=1;
maxItr=10^3;


parfor j=1:100
    j
[Ms2,Mm2,mean12]=MCCovHetro( sigt, 1, box_min,box_max,box_bin, l, v,1,1,maxItr,lambda,doCBS,2,2,ampfunc);
tMs(:,:,:,:,j)=Ms2; tMm(:,:,:,:,j)=Mm2; tmean1(:,:,j)=mean12;
end
Ms2=mean(tMs,5);
Mm2=mean(tMm,5);
mean12=mean(tmean1,3);
M2=Ms2+Mm2;



parfor j=1:100
    j
[Ms,Mm,mean1]=MCCovHetro( sigt, 1, box_min,box_max,box_bin, l, v,1,1,maxItr,lambda,doCBS,smpFlg,2,ampfunc);
tMs(:,:,:,:,j)=Ms; tMm(:,:,:,:,j)=Mm; tmean1(:,:,j)=mean1;
end
Ms=mean(tMs,5);
Mm=mean(tMm,5);
mean1=mean(tmean1,3);
M=Ms+Mm;



maxItrSi=10^2;
maxItrS=10^4;
uL=zeros(Nv,Nl,maxItrS);
if 1
parfor j=1:maxItrS
    j
    uL(:,:,j)=MCsampleHetro( sigt, 1, box_min,box_max,box_bin, l, v,1,1,maxItrSi,lambda,doCBS,2,2,ampfunc);
    
end
end
for j1=1:Nl
    for j2=1:Nl
       MU(:,:,j1,j2)=1/maxItrS*squeeze(uL(:,j1,:))*squeeze(uL(:,j2,:))';
    end
end
theta_v_d=rad2deg(theta_v);
theta_l_d=rad2deg(theta_l);

for j1=1%:Nl%:Nl
    for j2=j1:Nl
     
        ds=find(abs(theta_v_d-theta_l_d(j1))<0.001)-find(abs(theta_v_d-theta_l_d(j2))<0.001);
     
       figure, plot(abs(diag( M(:,:,j1,j2),ds))); hold on
       plot(abs(diag( M2(:,:,j1,j2),ds))); 
       plot(abs(diag( covU(:,:,j1,j2 ),ds)),'r');
          plot(abs(diag( MU(:,:,j1,j2),ds)));
        legend('MC','MC exp ','mu-diff','MC smp')
       title(sprintf('theta_1=%d, theta_2=%d',round(rad2deg(theta_l(j1))),round(rad2deg(theta_l(j2)))))
      
   
    end
end
figure, hold on, plot(theta_v_d,abs(meanU(:,1))),plot(theta_v_d,abs(mean1(:,1)))
figure, hold on, plot(theta_v_d,abs(meanU(:,3))),plot(theta_v_d,abs(mean1(:,3)))
