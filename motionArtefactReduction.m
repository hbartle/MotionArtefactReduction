%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Motion Artefact Reduction Algorithm
%
%   Wearable Electronic Devices Course E17
%   M.Sc. Electrical Engineering/Computer Engineering
%   Aarhus University
%   2017
%
%   Tongtong Jiang, Petr Kryze, Hannes Bartle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Measurement Data

% Delta-Sigma-ADC (Electrode)
load('data/steady/delsig_test1.mat');
% SAR1 (Microphone 1)
load('data/steady/sar1_test1.mat');
% SAR2 (Microphone 2)
load('data/steady/sar2_test1.mat');

%% Parameter Definition

% NLMS Filter Order
filterOrder = 10;
% NLMS Step Size
stepSize = 0.1;

% Power Estimation Weight Factor
beta = 0.3;
% Weight Adaption Sensitivity
gamma = 0.1;

%% Signal Initialization


%% Processing Loop
for i = 1:length(delsig)
   
    % High Pass Filter to remove offset
    xElectrode(i) = highPassFilter(delsig(i));
    xMic1(i) = highPassFilter(sar1(i));
    xMic2(i) = highPassFilter(sar2(i));
    
    % Adaptive FIR
    uMic1(i) = nlmsCoefficients1 * xMic1(i-filterOrder:i);
    uMic2(i) = nlmsCoefficients2 * xMic2(i-filterOrder:i);
    
    % Construct the error signal
    error1(i) = xElectrode(i) - uMic1(i);
    error2(i) = xElectrode(i) - uMic2(i);
    
    % Update the Adaptive FIR Filter Coefficients
    nlmsCoefficients1 = updateNLMS(nlmsCoefficients1,xMic1(i-filterOrder:i),...
                                   error1(i),stepSize);
    nlmsCoefficients2 = updateNLMS(nlmsCoefficients2,xMic2(i-filterOrder:i),...
                                   error2(i),stepSize);

    % Estimate Power of the Prediction Error
    powerError1(i) = predictionErrorEstimation(powerError1(i-1),...
                                               error1(i), beta);
    powerError2(i) = predictionErrorEstimation(powerError2(i-1),...
                                               error2(i), beta);
    
    % Update the Signal Combination Weight
    alpha(i) = weightAdaption(alpha(i-1),powerError1(i),...
                              powerError2(i),gamma);
    
    % Calculate the Output Signal
    output(i) = alpha(i) * error1(i) + (1-alpha(i)) * error2(i);
end


%% Post Processing and Plotting
inputSignalFigure = figure('NumberTitle','off','Name','Input Signals');
subplot(3,1,1);
plot(delsig);
subplot(3,1,2);
plot(sar1);
subplot(3,1,3);
plot(sar2);



