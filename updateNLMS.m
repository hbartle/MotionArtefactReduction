function [ newCoefficients ] = updateNLMS( oldCoefficients, inputSignal,...
                                           error, stepSize )
%updateNLMS Update the filter coefficients of the normalized
%least-mean-square filter


DEN = inputSignal*inputSignal' + 0.01;

newCoefficients = oldCoefficients + stepSize * inputSignal * error / DEN;

end

