
%% Generate ensamble

B = 199;
[ensamble, CovMat, L, mean] = genRealizations(B);
plotEnsamble(ensamble);

%% Forecast/Predict travel time

load('travelTimeData');
ensambleTravelTimes = forecastTravelTime(ensamble, B);

%% Calc empirical covariances

travelTimesMeans = -ones(50, 1);
travelTimesVariances = -ones(50, 1);

ensambleMeans = -ones(100, 1);
% ensambleVariances = -ones(100, 1);

for i = 1:50
    travelTimesMeans(i) = sum(ensambleTravelTimes(i, :))/B;
end

for i = 1:100
    ensambleMeans(i) = sum(ensamble(i, :))/B;
end

for i = 1:50
    travelTimesVariances(i) = sum((ensambleTravelTimes(i, :) - travelTimesMeans(i)).^2)/(B - 1);
end

% for i = 1:100
%     ensambleVariances(i) = sum((ensamble(i, :) - ensambleMeans(i)).^2)/(B - 1);
% end

%%

figure(2); hold off;
for i = 1:B
    plot(1:50, ensambleTravelTimes(:,i)); hold on;
end
plot(1:50, travelTimesMeans, 'b', 'LineWidth', 5);
plot(1:50, travelTimesMeans + travelTimesVariances, 'b', 'LineWidth', 5);
plot(1:50, travelTimesMeans - travelTimesVariances, 'b', 'LineWidth', 5);
grid on;

%%

% Assimilate data
