% Calculating dB level of # turbines at
% d_i distance for v_i wind speed. 

Noise(21,1.8)

function[NoisePerBuilding] = Noise(R,omega)
Dist = importdata("Distance.txt");

R = R;          % Radius [m]
omega = omega;

%Vtip = omega*R;      % Tip speed, m/s

%A = [           % A = length of # turbines - [distance, rotational speed]
%    [250,omega]
%    [160,omega]
%    [165,omega]
%    [230,omega]
%    [320,omega]
%    [410,omega]
%    [555,omega]
%    [515,omega]
%    [495,omega]
%    [510,omega]
%    ];
Itotal = 0;
alpha = 0.005;   % Frequency-dependent sound absorption coefficient, broadband estimate [dB(A)/m]
NoisePerBuilding = [];

for j = 1:length(Dist(1,:))
    A = [Dist(:,j), omega*ones(10,1)];
for i = 1:length(A(:,1))                            % For all of turbines in A
    Vtip = A(i,2)*R;
    Lwa = 50*(log10(Vtip))+10*(log10(2*R))-4;     % Calculate power at source            [W]
    Lp = Lwa - 10*log10(2*pi*(A(i,1))^2)- alpha*A(i,1);  % Calculate sound pressure level at distance A(i,1)
    I = LevelToIntensity(Lp);
    Itotal = Itotal + I;                        % Sum all intensities                       [W/m^2]

    NoiseindB = IntensityToLevel(Itotal);
    
end
NoisePerBuilding = [NoisePerBuilding, NoiseindB];
%S = 'Sound level in dB:';
%disp(S)
%disp(j)
%disp(NoiseindB)
end

% Transform intensity to dB
% and print value
%NoiseindB = IntensityToLevel(Itotal);
%S = 'Sound level in dB:';
%disp(S)
%disp(NoiseindB)

% Transform dB level to intensity
function L2I = LevelToIntensity(NoiseLevelIndB)
    I0 = 10.^(-12);                     % This is the treshold hearing intensity, matching 0 dB
    L2I = I0*10^(NoiseLevelIndB/10);    % I0 is reference value
end

% Transform intensity to dB level
function I2L = IntensityToLevel(Intensity)
    I0 = 10.^(-12);                 % This is the treshold hearing intensity, matching 0 dB
    I2L = 10*log10(Intensity/I0);   % Transform to log scale
end
end