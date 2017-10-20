function [ newAlpha ] = weightAdaption( oldAlpha, estimatedPower1,...
                                        estimatedPower2, gamma )
%WEIGHTADAPTION Calculate the weight of each predicted signal using the
%prediction error estimation.

newAlpha =oldAlpha - gamma * sign(estimatedPower1 -estimatedPower2);
if newAlpha > 1
    newAlpha = 1;
elseif newAlpha < 0
    newAlpha = 0;
end
end

