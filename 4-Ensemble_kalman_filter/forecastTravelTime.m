function ensambleTravelTimes  = forecastTravelTime(ensamble, B)

ensambleTravelTimes = -ones(50, B);
tau = 0.1;

% iterate throung ensambles
for i = 1:B
    
    
    % iterate through receivers
    for j = 1:50
    
        % add up all slownesses down to depth of receiver j, of ensamble i
        sumSlowness = 0;
        for k = 1:(50 + j)
            sumSlowness = sumSlowness + ensamble(k, i);
        end
        
        angleToReceiver = atan2(40, 50 + j);
        noise = tau*randn(1);
        
        % travel time to receiver j, of ensamble i
        ensambleTravelTimes(j, i) = sumSlowness/angleToReceiver + noise;
    end    
    
end

end
