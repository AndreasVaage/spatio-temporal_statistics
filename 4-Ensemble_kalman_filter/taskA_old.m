%% Generate ensamble

B = 199;
layers = 100;
load('travelTimeData');

[ensamble, CovMat, L, MU] = genRealizations(B, layers);
plotEnsamble(2, ensamble);

%% Forecast/Predict travel time


ensambleTravelTimes = forecastTravelTime(ensamble, B, 1:50);

%% Calc empirical covariances

ensambleTravelTimeCovariances = -ones(layers, 50);

travelTimesMeans = mean(ensambleTravelTimes,2);
travelTimesVariances = var(ensambleTravelTimes,1,2);

ensambleMeans = mean(ensamble,2);

for j = 1:50 % iterating throung sensors
    for k = 1:layers 
         C = cov(ensamble(k,:),ensambleTravelTimes(j,:));                  % Covariance matrix
         ensambleTravelTimeCovariances(k,j) = C(1,2);
    end
end

%%

figure(1); hold off;
for b = 1:B
    plot(1:50, ensambleTravelTimes(:,b)); hold on;
end
plot(1:50, travelTimesMeans, 'b', 'LineWidth', 5);
plot(1:50, travelTimesMeans + sqrt(travelTimesVariances), 'b', 'LineWidth', 5);
plot(1:50, travelTimesMeans - sqrt(travelTimesVariances), 'b', 'LineWidth', 5);
grid on;

%% Assimilate data

ensamble_assimilate = ensamble;

% Assimilate sensor j = 1
j = 50;
K = ensambleTravelTimeCovariances(:, j) / travelTimesVariances(j);       % Kalman gain for sensor j and realization i
h =  findobj('type','figure');
n_f = length(h)+1;
figure(n_f);
hold on
for b = 1:B
    ensamble_assimilate(:, b) = ensamble(:, b) + K *(travelTimeData(j) - ensambleTravelTimes(j,b));
end
plotEnsamble(n_f, ensamble_assimilate);

% Assimilate the other sensors 
for j = 49:-1:1
    K = ensambleTravelTimeCovariances(:, j) / travelTimesVariances(j);       % Kalman gain for sensor j and realization i
    for b = 1:B
        ensamble_assimilate(:, b) = ensamble_assimilate(:, b) + K *(travelTimeData(j) - ensambleTravelTimes(j,b));
    end
    
    plotEnsamble(n_f, ensamble_assimilate);
    pause(0.2)
end

plotEnsamble(j, ensamble_assimilate);

