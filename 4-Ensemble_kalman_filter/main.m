clear; close all
%% Generate ensamble

B = 199;
[ensamble, CovMat, L, MU] = genRealizations(B);
plotEnsamble(2, ensamble);

%% Forecast/Predict travel time

load('travelTimeData');
ensambleTravelTimes = forecastTravelTime(ensamble, B);

%% Calc empirical covariances

ensambleTravelTimeCovariances = -ones(100, 50);

travelTimesMeans = mean(ensambleTravelTimes,2);
travelTimesVariances = var(ensambleTravelTimes,1,2);

ensambleMeans = mean(ensamble,2);

for j = 1:50 % iterating throung sensors
    for k = 1:100 
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
K = ensambleTravelTimeCovariances(:, j)/travelTimesVariances(j);       % Kalman gain for sensor j and realization i
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
    K = ensambleTravelTimeCovariances(:, j)/travelTimesVariances(j);       % Kalman gain for sensor j and realization i
    for b = 1:B
        ensamble_assimilate(:, b) = ensamble_assimilate(:, b) + K *(travelTimeData(j) - ensambleTravelTimes(j,b));
    end
    
    plotEnsamble(n_f, ensamble_assimilate);
    pause(0.2)
end
%%
figure(n_f+1)
[yMean, yCIpercen] = CredInt(ensamble_assimilate',0.9);
plot(yMean, 1:100,'k'); hold on;
plot(yMean+yCIpercen,1:100,'--k')
ax = gca;
ax.YDir = 'reverse';
