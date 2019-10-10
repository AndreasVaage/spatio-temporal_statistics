clear; close all

%% A- Ensamble Kalman Filter Correct solution

B = 199;
depth = 100;
depths = 1:depth;
load('travelTimeData');

[slownessEnsamble, CovMat, L, MU] = genRealizations(B, depth);
plotEnsamble(1, slownessEnsamble);

% preallocate
slownessTravelTimeCovariance = -ones(depth,1);
slownessEnsamble_assimilate = slownessEnsamble;

for j = 1:50 % iterate through sensors
    
    % forecast travel time for sensor j; 
        % B different values, one for each ensamble
    sensorTravelTimeEnsamble = forecastTravelTime(slownessEnsamble_assimilate, B, j);
    
    % calculate variance of travel time for sensor j
    travelTimesVariance = var(sensorTravelTimeEnsamble);
    
    % calculate covariance of travel and with the slowness ensambles
    for k = depths
         C = cov(slownessEnsamble_assimilate(k,:), sensorTravelTimeEnsamble); % 2x2 Covariance matrix
         slownessTravelTimeCovariance(k) = C(1,2);
    end
    
    % calc kalman gain and assimilate data into ensamble
    K = slownessTravelTimeCovariance / travelTimesVariance;

    for b = 1:B
        slownessEnsamble_assimilate(:,b) = slownessEnsamble_assimilate(:,b) + K *(travelTimeData(j) - sensorTravelTimeEnsamble(b));
    end
    
    if j == 25
        slownessEnsamble_j_25 = slownessEnsamble_assimilate;
    end

    % plotEnsamble(2, slownessEnsamble_assimilate);
    % pause(0.1);
    
end

%%

slownessMean_j_25 = mean(slownessEnsamble_j_25, 2);
slownessVariance_j_25 = var(slownessEnsamble_j_25, 1, 2);

slownessMean = mean(slownessEnsamble_assimilate, 2);
slownessVariance = var(slownessEnsamble_assimilate, 1, 2);

std_90_j_25 = sqrt(slownessVariance_j_25)*norminv(0.9);
CI_j_25 = [slownessMean_j_25 - std_90_j_25, slownessMean_j_25 + std_90_j_25];

std_90 = sqrt(slownessVariance)*norminv(0.9);
CI = [slownessMean - std_90, slownessMean + std_90];

figure(3); hold off;
orange = [1, 0.5, 0]; blue = [0, 0.5, 1];
a1 = plot(slownessMean_j_25, depths, 'Color', blue, 'DisplayName', 'Estimated Slowness after 25 samples'); hold on; grid on;
a2 = plot(CI_j_25, 1:100, '--', 'Color', blue, 'DisplayName', '$90 \%$ CI');

b1 = plot(slownessMean, depths, 'Color', orange, 'DisplayName', 'Estimated Slowness all Sensors'); hold on; grid on;
b2 = plot(CI, 1:100, '--', 'Color', orange, 'DisplayName', '$90 \%$ CI');

ax = gca;
ax.YDir = 'reverse';
legend([a1, a2(1), b1, b2(1)], 'Location', 'best', 'interpreter', 'latex', 'FontSize', 10);
xlabel('Slowness [ms/m]', 'interpreter', 'latex', 'FontSize', 15);
ylabel('Depth [m]', 'interpreter', 'latex', 'FontSize', 15);



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

figure(4);
std_80_percent = (sqrt(diag(Sigma)))*norminv(0.8);
CI = [mu - std_80_percent, mu + std_80_percent];
plot(mu, 1:100,'k'); hold on;
plot(CI,1:100,'--k')
ax = gca;
ax.YDir = 'reverse';
legend("Estimated Slowness"," $80 \%$ CI",'Location','northwest','interpreter', 'latex', 'FontSize', 15);
xlabel("Slowness [ms/m]",'interpreter', 'latex', 'FontSize', 15);
ylabel("Depth [m]",'interpreter', 'latex', 'FontSize', 15);

