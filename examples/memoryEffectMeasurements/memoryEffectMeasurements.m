%% Define the expirements we want to run

% Delta between light to view (in order to avoid CBS)
theta0 = 0.002;

% light
theta = deg2rad((0:0.001:0.02) + theta0);

% phantom width
L = 2000;

% HG parameters
g = 0;
g0 = 0;

% OD
OD = 10;

% backward scattering
bws = false;

% view which will be squee
viewsVec = deg2rad(-1:0.1:-0.1);
vecNum = numel(viewsVec);

% sectors of cov matrix
sec_1 = 1:1:vecNum;
sec_2 = vecNum+1:1:2*vecNum;

% results
c_measured = zeros(1,numel(theta));

%% Run the expirements
    
% prepare the config
boxTargetArea = boxArea( ...
    1 ,             ... wavelength
    L/OD,           ... MFP
    [0,L],          ... z
    [-50000,50000], ... half of x axis
    [-50000,50000]  ... half of x axis
);

scatter = HGScatter(g,g0);

tic
% Run seperatly for each two directions
for idx1 = 1:1:numel(theta)
    theta_1 = theta(1);
    theta_2 = theta(idx1);

    dtheta = theta_2 - theta_1;
    phi_1 = viewsVec;
    phi_2 = viewsVec - dtheta;

    % If backward scattering
    if(bws)
        phi_1 = phi_1 + pi;
        phi_2 = phi_2 + pi;
    end
        
    % Run expirements
    mulres = scmc(boxTargetArea,                  ...
        farFieldSource([phi_1,phi_2],0),     ...
        farFieldSource([theta_1,theta_2],0), ...
        scatter,1e3, ... 1e5 for more accurate result
        'CBS',false, ... The light and the view are never the same direction
        'uniformFirstScatter',false, ... for more effitient sampling
        'parforIters', 12);

     % Get the intensity correlation
     sigma = zeros(2);
     sigma(1,1) = mean(diag(mulres.C(sec_1,sec_1,1,1)));
     sigma(2,2) = mean(diag(mulres.C(sec_2,sec_2,2,2)));
     
     % Check wich diagonal has higher values
     diagA = mean(diag(mulres.C(sec_2,sec_1,1,2)));
     diagB = mean(diag(mulres.C(sec_1,sec_2,1,2)));
     
     if(abs(diagA) > abs(diagB))
        sigma(1,2) = mean(diag(mulres.C(sec_2,sec_1,1,2)));
        sigma(2,1) = mean(diag(mulres.C(sec_1,sec_2,2,1)));
     else
        sigma(1,2) = mean(diag(mulres.C(sec_1,sec_2,1,2)));
        sigma(2,1) = mean(diag(mulres.C(sec_2,sec_1,2,1)));
     end

     sigma_I = zeros(2);
     sigma_I(1,1) = (sigma(1,1))^2;
     sigma_I(2,2) = (sigma(2,2))^2;
     sigma_I(1,2) = sigma(1,2) * sigma(2,1);

     sigma_I = real(sigma_I);

     % Calculate the memory effect
     c_measured(idx1) = (sigma_I(1,2)/(sqrt(sigma_I(2,2)*sigma_I(1,1))));
end
toc

%% Compare between the calculated and theoreric ME graphs
% The theoretic ME is true for the diffusive region
dtheta = theta - theta(1);
c_theory = (2*pi*L*dtheta./sinh(2*pi*L*dtheta)).^2;
c_theory(dtheta == 0) = 1;

figure
hold on
plot(rad2deg(dtheta),c_measured,'lineWidth',2);
plot(rad2deg(dtheta),c_theory,'k--','lineWidth',2);

xlabel('\theta[deg]');
ylim([0,1]);
legend('Measured','Theoretic')