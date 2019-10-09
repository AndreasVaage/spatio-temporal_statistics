
%% Generate ensamble

B = 199;
[ensamble, CovMat, L, mean] = genRealizations(B);
plotEnsamble(2, ensamble);

%% Forecast/Predict travel time

load('travelTimeData');
ensambleTravelTimes = forecastTravelTime(ensamble, B);

%% Calc empirical covariances

travelTimesMeans = -ones(50, 1);
travelTimesVariances = -ones(50, 1);

ensambleMeans = -ones(100, 1);
ensambleTravelTimeCovariances = -ones(100, 50);

for i = 1:50
    travelTimesMeans(i) = sum(ensambleTravelTimes(i, :))/B;
end
for i = 1:50
    travelTimesVariances(i) = sum((ensambleTravelTimes(i, :) - travelTimesMeans(i)).^2)/(B - 1);
end


for i = 1:100
    ensambleMeans(i) = sum(ensamble(i, :))/B;
end


for k = 1:100 % iterating through layers
    
    for j = 1:50 % iterating throung sensors
        ensambleTravelTimeCovariances(k, j) = sum((ensamble(k, :) - ensambleMeans(k)).*(ensambleTravelTimes(j, :) - travelTimesMeans(j)))/(B - 1);
    end
    %ensambleVariances(i) = sum((ensamble(i, :) - ensambleMeans(i)).^2)/(B - 1);
end

%%

figure(1); hold off;
for i = 1:B
    plot(1:50, ensambleTravelTimes(:,i)); hold on;
end
plot(1:50, travelTimesMeans, 'b', 'LineWidth', 5);
plot(1:50, travelTimesMeans + travelTimesVariances, 'b', 'LineWidth', 5);
plot(1:50, travelTimesMeans - travelTimesVariances, 'b', 'LineWidth', 5);
grid on;

%% Assimilate data

ensamble_assimilate = ensamble;

% Assimilate sensor j = 1
j = 1;
for i = 1:B
    ensamble_assimilate(:, i) = ensamble(:, i) + ensambleTravelTimeCovariances(:, j)/travelTimesVariances(j)*(travelTimeData(j) - ensambleTravelTimes(j,i));
end
plotEnsamble(3, ensamble_assimilate);

% Assimilate the other sensors j = 1
for j = 2:50
    for i = 1:B
        ensamble_assimilate(:, i) = ensamble_assimilate(:, i) + ensambleTravelTimeCovariances(:, j)/travelTimesVariances(j)*(travelTimeData(j) - ensambleTravelTimes(j,i));
    end
end
plotEnsamble(j, ensamble_assimilate);

