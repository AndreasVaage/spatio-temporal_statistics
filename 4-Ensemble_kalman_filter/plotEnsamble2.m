function [] = plotEnsamble2(j, slownessEnsamble)

depths = 1:100;

slownessMean = mean(slownessEnsamble, 2);
slownessVariance = var(slownessEnsamble, 1, 2);

std_90 = sqrt(slownessVariance)*norminv(0.9);
CI = [slownessMean - std_90, slownessMean + std_90];

[yMean, yCIpercen] = CredInt(slownessEnsamble',0.95);
inBetween = [(yMean + yCIpercen(1,:)), fliplr((yMean + yCIpercen(2,:)))];

plotEnsamble(j, slownessEnsamble);
p1 = plot(yMean, depths, 'k', 'linewidth', 1, 'DisplayName', 'Mean ensabled slowness'); hold on;
p2 = plot(CI, depths, '--k', 'linewidth', 1, 'DisplayName', '$10\%$ and $90\%$ uncertainty bounds');
fplt = fill(inBetween, [depths, fliplr(depths)], 'b');
set(fplt, 'facealpha', .2, 'DisplayName', '$95\%$ confidence of mean');
ax = gca;
ax.YDir = 'reverse';
legend([p1, fplt, p2(1)], 'location', 'northwest', 'interpreter', 'latex', 'FontSize', 10);
end

