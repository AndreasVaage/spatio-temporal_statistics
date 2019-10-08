
%% Generate ensamble

B = 200;
[ensamble, CovMat, L, mean] = genRealizations(B);
plotEnsamble(ensamble);

%% Forecast/Predict travel time

load('travelTimeData');
ensambleTravelTimes = forecastTravelTime(ensamble, B);

% Calc empirical covariances
% Assimilate data
