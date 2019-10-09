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

plotEnsamble(j, ensamble_assimilate);

%% C - Kalman filter

Sigma = CovMat;
mu = MU;
tau = 0.1;
g = zeros(1,100);
for j = 1:50
    angleToReceiver = atan2(40, 50 + j);
    g(1:50+j) = 1/cos(angleToReceiver);
    g(50+j:100) = 0;
    
    K = Sigma*g'/(g*Sigma*g' + tau^2);
    mu = mu + K*(travelTimeData(j) - g*mu);
    Sigma = Sigma - K*g*Sigma;
end

figure
std_80_percent = (sqrt(diag(Sigma)))*norminv(0.8);
CI = [mu - std_80_percent, mu + std_80_percent];
plot(mu, 1:100,'k'); hold on;
plot(CI,1:100,'--k')
ax = gca;
ax.YDir = 'reverse';
legend("Estimated Slowness"," $80 \%$ CI",'Location','northwest','interpreter', 'latex', 'FontSize', 15)
xlabel("Slowness [ms/m]",'interpreter', 'latex', 'FontSize', 15)
ylabel("Depth [m]",'interpreter', 'latex', 'FontSize', 15)

