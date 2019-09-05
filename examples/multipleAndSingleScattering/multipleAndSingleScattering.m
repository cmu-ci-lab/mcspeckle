%% Check dependence of single vs. multiple scattering
% Measure the correlation in multiple scattering and single scattering
% cases, where the light and view are rotated around the z axis

span = 22; % how l1 source is far from z axis (in degrees)
alpha = 0:3:90; % the rotation directions (in degrees)

% Build the target area
boxTargetArea = boxArea( ...
    1 ,       ... wavelength
    20,       ... MFP
    [-50,50], ... z
    [-50,50], ... x
    [-50,50]  ... y
);

scatter = isotropicScatter;

singleCorr = zeros(1,numel(alpha));
mulCorr = zeros(1,numel(alpha));
totalCorr = zeros(1,numel(alpha));

%% Measure correlation

l1 = [-sind(span), 0,  cosd(span)];
v1 = [-sind(span), 0, -cosd(span)]; % minus for backward scattering

tic
for a = 1:1:numel(alpha)
    % Build the rotation matrix around z axis
    R_alpha = [ ...
        cosd(alpha(a)), -sind(alpha(a)), 0; ...
        sind(alpha(a)),  cosd(alpha(a)), 0; ...
        0                0               1];
    
    l2=(R_alpha*l1')';
    v2=(R_alpha*v1')';
    
    % Solve
    mulres = scmc(boxTargetArea,   ...
        farFieldSource([v1',v2']), ...
        farFieldSource([l1',l2']), ...
        scatter,1e3, ...
        'mean', false, ...
        'singleScattering', true, ...
        'parforIters', 12);
    
    singleCorrMatrix = mulres.Csingle(:,:,1,2);
    totalCorrMatrix = mulres.C(:,:,1,2);
    [~,maxCorrIdx] = max(abs(totalCorrMatrix(:)));
    
    singleCorr(a) = abs(singleCorrMatrix(maxCorrIdx));
    mulCorr(a) = abs(totalCorrMatrix(maxCorrIdx) - singleCorrMatrix(maxCorrIdx));
    totalCorr(a) = abs(totalCorrMatrix(maxCorrIdx));
    
end
toc

maxCorr = max(totalCorr);

%% Comapare between single scattering the multiple scattering
figure;
hold on
plot(alpha,singleCorr/maxCorr,'lineWidth',2);
plot(alpha,mulCorr/maxCorr,'lineWidth',2);
plot(alpha,totalCorr/maxCorr,'lineWidth',2);

legend('Single Scattering','Multiple Scattering','Total Scattering')
xlabel('Rotation [deg]');
ylabel('Correlation');