function [ newPower] = predictionErrorEstimation(oldPower, signal, beta )
%PREDICTIONERRORESTIMATION Estimate the error in the predicted (filtered)
%signals
newPower = (1-beta)*oldPower + beta*signal.^2;

end
