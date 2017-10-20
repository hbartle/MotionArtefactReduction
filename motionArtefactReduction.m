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
filterOrderNLMS = 30;
% NLMS Step Size
stepSize = 0.001;

% HP Filter Order
filterOrderHP = 6;
% HP Passpand Frequency;
passBandFrequencyHP = 2.8;

% Notch Filter Center Frequency
notchFrequency = 50/(sampleRate/2);
% Notch Filter Bandwidth
notchBandwidth = notchFrequency/30;


% Power Estimation Weight Factor
beta = 0.8;
% Weight Adaption Sensitivity
gamma = 0.3;

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
                    'PassbandFrequency',passBandFrequencyHP,'PassbandRipple',1e-6, ...
                    'SampleRate',sampleRate);

nlmsCoefficients1 = zeros(1,filterOrderNLMS);                         
nlmsCoefficients2 = zeros(1,filterOrderNLMS);  

[notchFiltB, notchFiltA] = iirnotch(notchFrequency,notchBandwidth);

%% Processing Loop
for i = 1:length(delsig)
   
    % High Pass Filter to remove offset
    hpElectrode(1:i) = filter(hpFilt,delsig(1:i));
    hpMic1(1:i) = filter(hpFilt,sar1(1:i));
    hpMic2(1:i) = filter(hpFilt,sar2(1:i));
%     hpElectrode(1:i) = delsig(1:i) - mean(delsig(1:i));
%     hpMic1(1:i) = sar1(1:i)- mean(sar1(1:i));
%     hpMic2(1:i) = sar2(1:i)- mean(sar2(1:i));
    
    % Notch Filter to remove 50Hz 
    notchElectrode(1:i) = filter(notchFiltB,notchFiltA,hpElectrode(1:i));
    notchMic1(1:i) = filter(notchFiltB,notchFiltA,hpMic1(1:i));
    notchMic2(1:i) = filter(notchFiltB,notchFiltA,hpMic2(1:i));
    
    % Adaptive FIR
    uMic1(1:i) = filter(nlmsCoefficients1,1,notchMic1(1:i));
    uMic2(1:i) = filter(nlmsCoefficients2,1,notchMic2(1:i));
    
    % Construct the error signal
    error1(i) = notchElectrode(i) - uMic1(i);
    error2(i) = notchElectrode(i) - uMic2(i);
    
    
    % Update the Adaptive FIR Filter Coefficients
    if i < filterOrderNLMS 
        % Adapt the Filter Input Vector
        filterInput1 = [hpMic1(1:i), zeros(1,filterOrderNLMS-i)];
        filterInput2 = [hpMic2(1:i), zeros(1,filterOrderNLMS-i)];

        nlmsCoefficients1 = updateNLMS(nlmsCoefficients1,...
                                       filterInput1,...
                                       error1(i),stepSize);
        nlmsCoefficients2 = updateNLMS(nlmsCoefficients2,...
                                       filterInput2,...
                                       error2(i),stepSize);
    else
        nlmsCoefficients1 = updateNLMS(nlmsCoefficients1,...
                                       hpMic1(i-filterOrderNLMS+1:i),...
                                       error1(i),stepSize);
        nlmsCoefficients2 = updateNLMS(nlmsCoefficients2,...
                                       hpMic2(i-filterOrderNLMS+1:i),...
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
close all
plotMin = 900;
plotMax = 1900;
sampleTime = 1/sampleRate;
plotRange = plotMin*sampleTime:sampleTime:plotMax*sampleTime;

legendFontSize = 16;

%%%%%%%%%%%%%%%%%%%%%%
rawElectrodeFigure = figure('NumberTitle','off',...
                            'Name','Raw Electrode Signal',...
                            'units','normalized','outerposition',[0 0 1 1]);
plot(plotRange,delsig(plotMin:plotMax));
grid on
xlabel('Time [s]')
title('Raw Electrode Signal')
%%%%%%%%%%%%%%%%%%%%%%
rawMicFigure = figure('NumberTitle','off',...
                      'Name','Raw Microphone Signals',...
                      'units','normalized','outerposition',[0 0 1 1]);
plot(plotRange,sar1(plotMin:plotMax),...
     plotRange,sar2(plotMin:plotMax));
grid on
title('Microphones')
xlabel('Time [s]')
lgd=legend('Mic 1','Mic 2');
lgd.FontSize = legendFontSize;
%%%%%%%%%%%%%%%%%%%%%%

filtElectrodeFigure = figure('NumberTitle','off',...
                             'Name','Filtered Electrode Signals',...
                             'units','normalized','outerposition',[0 0 1 1]);
plot(plotRange,hpElectrode(plotMin:plotMax),...
     plotRange,notchElectrode(plotMin:plotMax));
grid on
xlabel('Time [s]')
title('Filtered Electrode Signal')
lgd=legend('HP','HP + Notch');
lgd.FontSize = legendFontSize;
%%%%%%%%%%%%%%%%%%%%%%

filtMicFigure = figure('NumberTitle','off',...
                       'Name','Filtered Microphone Signals',...
                       'units','normalized','outerposition',[0 0 1 1]);
plot(plotRange,uMic1(plotMin:plotMax),...
     plotRange,uMic2(plotMin:plotMax));
grid on
title('Filtered Microphone Signals')
xlabel('Time [s]')
lgd=legend('Mic 1','Mic 2');
lgd.FontSize = legendFontSize;
%%%%%%%%%%%%%%%%%%%%%%

alphaFigure = figure('NumberTitle','off',...
                     'Name','Weight Factor Alpha',...
                     'units','normalized','outerposition',[0 0 1 1]);
plot(plotRange,alpha(plotMin:plotMax));
grid on
title('Weight Factor (Alpha)')
xlabel('Time [s]')
%legend('Weight Factor (Alpha)')
%%%%%%%%%%%%%%%%%%%%%%

errorFigure = figure('NumberTitle','off',...
                     'Name','Error Signals',...
                     'units','normalized','outerposition',[0 0 1 1]);
plot(plotRange,error1(plotMin:plotMax),...
     plotRange,error2(plotMin:plotMax),... 
     plotRange,output(plotMin:plotMax));
grid on
title('Error Signals and Output')
xlabel('Time [s]')
lgd=legend('Error 1','Error 2','Final Signal');
lgd.FontSize = legendFontSize;
%%%%%%%%%%%%%%%%%%%%%%

outputFigure = figure('NumberTitle','off',...
                      'Name','Final Signal',...
                      'units','normalized','outerposition',[0 0 1 1]);plot(plotRange,hpElectrode(plotMin:plotMax),...    
plotRange,output(plotMin:plotMax));
grid on
title('Filter Output and Input  ')
xlabel('Time [s]')
lgd=legend('HP Electrode','Final Signal');
lgd.FontSize = legendFontSize;
%% Print Plots to EPS files
mkdir('results')
print(rawElectrodeFigure, strcat('results/rawElectrode_case_',int2str(dataCase)),'-depsc');
print(rawMicFigure, strcat('results/rawMic_case_',int2str(dataCase)),'-depsc');
print(filtElectrodeFigure, strcat('results/filtElectrode_case_',int2str(dataCase)),'-depsc');
print(filtMicFigure, strcat('results/filtMic_case_',int2str(dataCase)),'-depsc');
print(alphaFigure, strcat('results/alpha_case_',int2str(dataCase)),'-depsc');
print(errorFigure, strcat('results/error_case_',int2str(dataCase)),'-depsc');
print(outputFigure, strcat('results/output_case_',int2str(dataCase)),'-depsc');


