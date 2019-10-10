clear; close all

%% A - Ensamble Kalman Filter Correct solution

B = 200;
depth = 100;
depths = 1:depth;
load('travelTimeData');

[slownessEnsamble, CovMat, L, MU] = genRealizations(B, depth);

figure(1);
plotEnsamble(1, slownessEnsamble);

[slownessEnsamble_assimilate, slownessEnsamble_j_25] = enKF(slownessEnsamble, travelTimeData);

figure;
subplot(1,2,1);
plotEnsamble(1, slownessEnsamble_assimilate);
subplot(1,2,2);
plotEnsamble(2, X);

%%
slownessMean_j_25 = mean(slownessEnsamble_j_25, 2);
slownessVariance_j_25 = var(slownessEnsamble_j_25, 1, 2);

slownessMean = mean(slownessEnsamble_assimilate, 2);
slownessVariance = var(slownessEnsamble_assimilate, 1, 2);

std_90_j_25 = sqrt(slownessVariance_j_25)*norminv(0.9);
CI_j_25 = [slownessMean_j_25 - std_90_j_25, slownessMean_j_25 + std_90_j_25];

std_90 = sqrt(slownessVariance)*norminv(0.9);
CI = [slownessMean - std_90, slownessMean + std_90];

[yMean, yCIpercen] = CredInt(slownessEnsamble_assimilate',0.95);
inBetween = [(yMean+yCIpercen(1,:)), fliplr((yMean + yCIpercen(2,:)))];

figure(2);
subplot(1,2,1);
plotEnsamble(50, slownessEnsamble_assimilate);
p1 = plot(yMean, depths,'k','linewidth',2,'DisplayName','Mean ensabled slowness'); hold on;
p2 = plot(CI, depths, '--k','linewidth',2, 'DisplayName', '$10\%$ and $90\%$ uncertainty bounds');
fplt = fill(inBetween, [1:100, fliplr(1:100)], 'b');
set(fplt, 'facealpha', .2, 'DisplayName', '$95\%$ confidence of mean');
ax = gca;
ax.YDir = 'reverse';
legend([p1, fplt, p2(1)], 'location', 'northwest', 'interpreter', 'latex', 'FontSize', 12)

subplot(1,2,2);
orange = [1, 0.5, 0]; blue = [0, 0.5, 1];
a1 = plot(slownessMean_j_25, depths, 'Color', blue, 'DisplayName', 'Estimated Slowness after assimilating 25 sensors'); hold on; grid on;
a2 = plot(CI_j_25, 1:100, '--', 'Color', blue, 'DisplayName', '$10\%$ and $90\%$ uncertainty bound');
b1 = plot(slownessMean, depths, 'Color', orange, 'DisplayName', 'Estimated Slowness all Sensors'); hold on; grid on;
b2 = plot(CI, 1:100, '--', 'Color', orange, 'DisplayName', '$10\%$ and $90\%$ uncertainty bound');
ax = gca;
ax.YDir = 'reverse';
title('\textbf{Final and intermediate (j = 25) ensamble estimation}', 'interpreter', 'latex', 'FontSize', 18);
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
plot(mu, depths,'k'); hold on;
plot(CI, depths,'--k')
ax = gca;
ax.YDir = 'reverse';
legend("Estimated Slowness"," $10\%$ and $90\%$ uncertainty bounds",'Location','northwest','interpreter', 'latex', 'FontSize', 15);
xlabel("Slowness [ms/m]",'interpreter', 'latex', 'FontSize', 15);
ylabel("Depth [m]",'interpreter', 'latex', 'FontSize', 15);

