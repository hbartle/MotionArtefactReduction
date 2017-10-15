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
clear all
close all
clc 
%% Load Measurement Data

dataCase = 2;

if dataCase == 1
    % Delta-Sigma-ADC (Electrode)
    load('data/steady/delsig.mat');
    % SAR1 (Microphone 1)
    load('data/steady/sar1.mat');
    % SAR2 (Microphone 2)
    load('data/steady/sar2.mat');
    
elseif dataCase == 2
    % Delta-Sigma-ADC (Electrode)
    load('data/left_hand_moving/delsig.mat');
    % SAR1 (Microphone 1)
    load('data/left_hand_moving/sar1.mat');
    % SAR2 (Microphone 2)
    load('data/left_hand_moving/sar2.mat');
    
elseif dataCase == 3
    % Delta-Sigma-ADC (Electrode)
    load('data/both_hands_moving/delsig.mat');
    % SAR1 (Microphone 1)
    load('data/both_hands_moving/sar1.mat');
    % SAR2 (Microphone 2)
    load('data/both_hands_moving/sar2.mat');
end

% Sampling Frequency [Hz]
sampleRate = 183;

%% Parameter Definition

% NLMS Filter Order
filterOrderNLMS = 20;
% NLMS Step Size
stepSize = 0.001;

% HP Filter Order
filterOrderHP = 5;
% HP Passpand Frequency;
passBandFrequencyHP = 2;

% Power Estimation Weight Factor
beta = 0.4;
% Weight Adaption Sensitivity
gamma = 0.4;

%% Signal Initialization
powerError1 = 0;
powerError2 = 0;
alpha = 0.5;

% Fix Last Values
delsig = delsig(1:end-1);
sar1 = sar1(1:end-1);
sar2 = sar2(1:end-1);


%% Construct Filters
hpFilt = designfilt('highpassiir','FilterOrder',filterOrderHP, ...
                    'PassbandFrequency',passBandFrequencyHP,'PassbandRipple',0.1, ...
                    'SampleRate',sampleRate);

nlmsCoefficients1 = zeros(1,filterOrderNLMS);                         
nlmsCoefficients2 = zeros(1,filterOrderNLMS);                         
%% Processing Loop
for i = 1:length(delsig)
   
    % High Pass Filter to remove offset
    xElectrode(1:i) = filter(hpFilt,delsig(1:i));
    xMic1(1:i) = filter(hpFilt,sar1(1:i));
    xMic2(1:i) = filter(hpFilt,sar2(1:i));
    
    % Adaptive FIR
    uMic1(1:i) = filter(nlmsCoefficients1,1,xMic1(1:i));
    uMic2(1:i) = filter(nlmsCoefficients2,1,xMic2(1:i));
    
    % Construct the error signal
    error1(i) = xElectrode(i) - uMic1(i);
    error2(i) = xElectrode(i) - uMic2(i);
    
    
    % Update the Adaptive FIR Filter Coefficients
    if i < filterOrderNLMS 
        % Adapt the Filter Input Vector
        filterInput1 = [xMic1(1:i), zeros(1,filterOrderNLMS-i)];
        filterInput2 = [xMic2(1:i), zeros(1,filterOrderNLMS-i)];

        nlmsCoefficients1 = updateNLMS(nlmsCoefficients1,...
                                       filterInput1,...
                                       error1(i),stepSize);
        nlmsCoefficients2 = updateNLMS(nlmsCoefficients2,...
                                       filterInput2,...
                                       error2(i),stepSize);
    else
        nlmsCoefficients1 = updateNLMS(nlmsCoefficients1,...
                                       xMic1(i-filterOrderNLMS+1:i),...
                                       error1(i),stepSize);
        nlmsCoefficients2 = updateNLMS(nlmsCoefficients2,...
                                       xMic2(i-filterOrderNLMS+1:i),...
                                       error2(i),stepSize);
    end

    if i>1
        % Estimate Power of the Prediction Error
        powerError1(i) = predictionErrorEstimation(powerError1(i-1),...
            error1(i), beta);
        powerError2(i) = predictionErrorEstimation(powerError2(i-1),...
            error2(i), beta);
        
        % Update the Signal Combination Weight
        alpha(i) = weightAdaption(alpha(i-1),powerError1(i),...
            powerError2(i),gamma);
    end
        
    
    % Calculate the Output Signal
    output(i) = alpha(i) * error1(i) + (1-alpha(i)) * error2(i);
end


%% Post Processing and Plotting
plotMin = 500;
plotMax = 800;

inputSignalFigure = figure('NumberTitle','off','Name','Input Signals');
subplot(3,1,1);
plot(plotMin:plotMax,delsig(plotMin:plotMax));
title('Electrode Signal')
subplot(3,1,2);
plot(plotMin:plotMax,sar1(plotMin:plotMax));
title('Microphone 1')
subplot(3,1,3);
plot(plotMin:plotMax,sar2(plotMin:plotMax));
title('Microphone 2')

filteredSignalFigure = figure('NumberTitle','off','Name','Filtered Signals');
subplot(3,1,1);
plot(plotMin:plotMax,xMic1(plotMin:plotMax),...
     plotMin:plotMax,xMic2(plotMin:plotMax));
title('High Pass Filtered Microphone Signals')
legend('Mic 1','Mic 2')
subplot(3,1,2);
plot(plotMin:plotMax,uMic1(plotMin:plotMax),...
     plotMin:plotMax,uMic2(plotMin:plotMax));
title('NLMS Filtered Microphone Signals')
legend('Mic 1','Mic 2')

subplot(3,1,3);
plot(plotMin:plotMax,error1(plotMin:plotMax),...
     plotMin:plotMax,error2(plotMin:plotMax));
title('Error Signals')
legend('Error 1','Error 2')


autoArrangeFigures(1,2)

outputSignalFigure = figure('NumberTitle','off','Name','Output Signal');
subplot(2,1,1);
plot(plotMin:plotMax,alpha(plotMin:plotMax));
title('Weight Factor (Alpha)')
subplot(2,1,2);
plot(plotMin:plotMax,output(plotMin:plotMax));
title('Final Signal')



