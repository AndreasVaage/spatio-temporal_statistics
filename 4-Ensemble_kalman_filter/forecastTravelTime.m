function travelTimes  = forecastTravelTime(ensamble, B, sensors)

m = length(sensors);
travelTimes = -ones(m, B);
tau = 0.1;

% iterate throung ensambles
for i = 1:B
    
    % iterate through the given receivers/sensors
    for j = 1:m
    
        % add up all slownesses down to depth of receiver j, of ensamble i
        sumSlowness = sum(ensamble(1:50 + sensors(j),i));
        
        angleToReceiver = atan2(40, 50 + sensors(j));
        noise = tau*randn(1);
        
        % travel time to receiver j, of ensamble i
        travelTimes(j, i) = sumSlowness/cos(angleToReceiver) + noise;
    end    
    
end

end

