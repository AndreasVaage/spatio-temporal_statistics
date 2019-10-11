clear; close all

%% A - Ensamble Kalman Filter Correct solution

depth = 100;
depths = 1:depth;
load('travelTimeData');

[slownessEnsamble_200, CovMat, L, MU] = genRealizations(200);
[slownessEnsamble_200_j_50, slownessEnsamble_200_j_25] = enKF(slownessEnsamble_200, travelTimeData, 1);
[slownessEnsamble_200_j_50_rev, slownessEnsamble_200_j_25_rev] = enKF(slownessEnsamble_200, travelTimeData, 1);

slownessEnsamble_50 = genRealizations(50);
[slownessEnsamble_50_j_50, slownessEnsamble_50_j_25] = enKF(slownessEnsamble_50, travelTimeData, 1);
[slownessEnsamble_50_j_50_rev, slownessEnsamble_50_j_25_rev] = enKF(slownessEnsamble_50, travelTimeData, 1);

slownessEnsamble_400 = genRealizations(400);
[slownessEnsamble_400_j_50, slownessEnsamble_400_j_25] = enKF(slownessEnsamble_400, travelTimeData, 1);
[slownessEnsamble_400_j_50_rev, slownessEnsamble_400_j_25_rev] = enKF(slownessEnsamble_400, travelTimeData, 1);

%%
disp('200 Realizations');
figure(1);
subplot(1,3,1);
plotEnsamble2(1, slownessEnsamble_200);
subplot(1,3,2);
plotEnsamble2(25, slownessEnsamble_200_j_25);
subplot(1,3,3);
plotEnsamble2(50, slownessEnsamble_200_j_50);

disp('50 Realizations');
figure(2);
subplot(1,3,1);
plotEnsamble2(1, slownessEnsamble_50);
subplot(1,3,2);
plotEnsamble2(25, slownessEnsamble_50_j_25);
subplot(1,3,3);
plotEnsamble2(50, slownessEnsamble_50_j_50);

disp('400 Realizations');
figure(3);
subplot(1,3,1);
plotEnsamble2(1, slownessEnsamble_400);
subplot(1,3,2);
plotEnsamble2(25, slownessEnsamble_400_j_25);
subplot(1,3,3);
plotEnsamble2(50, slownessEnsamble_400_j_50);


disp('200 Realizations');
figure(4);
subplot(1,3,1);
plotEnsamble2(1, slownessEnsamble_200);
subplot(1,3,2);
plotEnsamble2(25, slownessEnsamble_200_j_25_rev);
subplot(1,3,3);
plotEnsamble2(50, slownessEnsamble_200_j_50_rev);

disp('50 Realizations');
figure(5);
subplot(1,3,1);
plotEnsamble2(1, slownessEnsamble_50);
subplot(1,3,2);
plotEnsamble2(25, slownessEnsamble_50_j_25_rev);
subplot(1,3,3);
plotEnsamble2(50, slownessEnsamble_50_j_50_rev);

disp('400 Realizations');
figure(6);
subplot(1,3,1);
plotEnsamble2(1, slownessEnsamble_400);
subplot(1,3,2);
plotEnsamble2(25, slownessEnsamble_400_j_25_rev);
subplot(1,3,3);
plotEnsamble2(50, slownessEnsamble_400_j_50_rev);




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

std_80_percent = (sqrt(diag(Sigma)))*norminv(0.8);
CI = [mu - std_80_percent, mu + std_80_percent];

% B = 200
slownessMean_200 = mean(slownessEnsamble_200_j_50, 2);
slownessVariance_200 = var(slownessEnsamble_200_j_50, 1, 2);
std_90_200 = sqrt(slownessVariance_200)*norminv(0.9);
CI_200 = [slownessMean_200 - std_90_200, slownessMean_200 + std_90_200];

% B = 50
slownessMean_50 = mean(slownessEnsamble_50_j_50, 2);
slownessVariance_50 = var(slownessEnsamble_50_j_50, 1, 2);
std_90_50 = sqrt(slownessVariance_50)*norminv(0.9);
CI_50 = [slownessMean_50 - std_90_50, slownessMean_50 + std_90_50];

% B = 400
slownessMean_400 = mean(slownessEnsamble_400_j_50, 2);
slownessVariance_400 = var(slownessEnsamble_400_j_50, 1, 2);
std_90_400 = sqrt(slownessVariance_400)*norminv(0.9);
CI_400 = [slownessMean_400 - std_90_400, slownessMean_400 + std_90_400];


figure(7); hold off;
orange = [1, 0.5, 0]; blue = [0, 0.5, 1];
other = [0.5, 0, 1];
k1 = plot(mu, depths,'k'); hold on;
k2 = plot(CI, depths,'--k');

a1 = plot(slownessMean_200, depths, 'Color', blue, 'DisplayName', 'Estimated Slowness 200'); hold on; grid on;
a2 = plot(CI_200, depths, '--', 'Color', blue, 'DisplayName', '$10\%$ and $90\%$ uncertainty bound');

b1 = plot(slownessMean_50, depths, 'Color', orange, 'DisplayName', 'Estimated Slowness 50'); hold on; grid on;
b2 = plot(CI_50, depths, '--', 'Color', orange, 'DisplayName', '$10\%$ and $90\%$ uncertainty bound');

c1 = plot(slownessMean_400, depths, 'Color', other, 'DisplayName', 'Estimated Slowness 400'); hold on; grid on;
c2 = plot(CI_400, depths, '--', 'Color', other, 'DisplayName', '$10\%$ and $90\%$ uncertainty bound');

ax = gca;
ax.YDir = 'reverse';
title('\textbf{Final and intermediate (j = 25) ensamble estimation}', 'interpreter', 'latex', 'FontSize', 18);
legend([k1, k2(1), a1, a2(1), b1, b2(1), c1, c2(1)], 'Location', 'best', 'interpreter', 'latex', 'FontSize', 10);


ax = gca;
ax.YDir = 'reverse';
legend("Estimated Slowness"," $10\%$ and $90\%$ uncertainty bounds",'Location','northwest','interpreter', 'latex', 'FontSize', 15);
xlabel("Slowness [ms/m]",'interpreter', 'latex', 'FontSize', 15);
ylabel("Depth [m]",'interpreter', 'latex', 'FontSize', 15);


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

function [ensamble, CovMat, L, mean_] = genRealizations(B)

layers = 100;

ensamble = -ones(layers, B);

CovMat = -ones(layers, layers);

sigma = 0.05;
mean_ = 0.5 - 0.001*(1:layers)';
eta = 0.1;

for i = 1:layers
    for j = 1:layers
        
        CovMat(i, j) = sigma^2 * (1 + eta*abs(i - j))*exp(-eta*abs(i - j));
        
    end
end

L = chol(CovMat)';

for i = 1:B
    ensamble(:, i) = mean_ + L*randn(layers, 1);
end


end


