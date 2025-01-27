function [] = plotEnsamble(j, ensamble)

s = size(ensamble);

if s(1) == 100 % one ensamble is a column
    A = s(2);
else % one ensamble is a row
    A = s(1);
end

hold off;
for i = 1:A
    plot(ensamble(:, i), 1:100, 'color', [0 0 0]+0.6); hold on;
end

ax = gca;
ax.YDir = 'reverse';
grid on;
title(['\textbf{Ensamble ' num2str(j) '}'], 'interpreter', 'latex', 'FontSize', 18);
xlabel('Slowness [ms/m]', 'interpreter', 'latex', 'FontSize', 15);
ylabel('Depth [m]', 'interpreter', 'latex', 'FontSize', 15);
xlim([0.3 0.6]);

end

