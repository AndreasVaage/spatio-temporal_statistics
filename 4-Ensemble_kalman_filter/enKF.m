function [slownessEnsamble_assimilate, slownessEnsamble_j_25] = enKF(slownessEnsamble, travelTimeData, order)
%ENKF function solves task (a)
%   given an ensamble and data, filters the ensamble using the data to
%   estimate parameters

s = size(slownessEnsamble);
depth = s(1);
B = s(2);

% preallocate
slownessTravelTimeCovariance = -ones(depth,1);
slownessEnsamble_assimilate = slownessEnsamble;

if order == 1
    sensors = 1:50;
else
    sensors = 50:-1:1;
end

for j = sensors % iterate through sensors
    
    % forecast travel time for sensor j; B times, one for each realization
    sensorTravelTimeEnsamble = forecastTravelTime(slownessEnsamble_assimilate, B, j);
    
    % calculate the travel time variance for sensor j
    travelTimesVariance = var(sensorTravelTimeEnsamble);
    
    % calculate covariance of travel with the slowness ensamble
    for k = 1:depth
         C = cov(slownessEnsamble_assimilate(k,:), sensorTravelTimeEnsamble); % 2x2 Covariance matrix
         slownessTravelTimeCovariance(k) = C(1,2);
    end
    
    % calc kalman gain and assimilate data into ensamble
    K = slownessTravelTimeCovariance / travelTimesVariance;
    for b = 1:B
        slownessEnsamble_assimilate(:,b) = slownessEnsamble_assimilate(:,b) + K*(travelTimeData(j) - sensorTravelTimeEnsamble(b));
    end
    
    if j == 25 % keeping track of the half-way/intermediate result
        slownessEnsamble_j_25 = slownessEnsamble_assimilate;
    end
    
end

end

