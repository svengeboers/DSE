% By Louis A. w/ input from Svennie

% Calculating dB level of # turbines at
% d_i distance for v_i wind speed. 

A = [           % A = length of # turbines - [distance,wind speed]
    [50,5.]
    [65,6]
    [60,4]
    ];
Itotal = 0;

% Calculate total local intensity (Itotal)
for i = 1:length(A(:,1))        % For all of turbines in A
    P = soundPower(A(i,2));     % Calculate power at source                 [W]
    I = P/(4*pi*(A(i,1))^2);    % Calculate intensity at distance A(i,1)    [W/m^2]   
    Itotal = Itotal + I;        % Sum all intensities                       [W/m^2]
end
   
% Transform intensity to dB
% and print value
NoiseindB = IntensityToLevel(Itotal);
S = 'Sound level in dB:';
disp(S)
disp(NoiseindB)

% Calculate sound power for certain wind speed v [W]
% From article of Sven/Louis A.
% Iu = 4*10^(-6.)*exp(0.2216*v) @ 0.4m 
% https://ijer.ut.ac.ir/article_1019_db4f5f4f140bb9dfa96c1aba59bcb60c.pdf
% Iu = 0.0034*exp(0.2533*v) @ 63m (?) - Vestas v52, 850 kW

function P = soundPower(v)
    % Iuncorrected = 0.0034*exp(0.2533*v);    % This data was measured at .4 m of source
    P = 0.0034*exp(0.2533*v);        % Power emitted
end

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