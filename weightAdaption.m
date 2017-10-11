function [ newAlpha ] = weightAdaption( oldAlpha, estimatedPower1,...
                                        estimatedPower2, gamma )
%WEIGHTADAPTION Calculate the weight of each predicted signal using the
%prediction error estimation.

newAlpha = oldAlpha - gamma * sgn(estimatedPower1 -estimatedPower2);

end

