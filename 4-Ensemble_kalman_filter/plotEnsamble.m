function [] = plotEnsamble(ensamble)


s = size(ensamble);

if s(1) == 100 % one ensamble is a column
    A = s(2);
else % one ensamble is a row
    A = s(1);
end
figure(1);
ax = gca;
ax.YDir = 'reverse';
grid on;
for i = 1:A
        plot(ensamble(:, i), 1:100); hold on;
end

title('\textbf{Ensamble-based prior}', 'interpreter', 'latex', 'FontSize', 18);
xlabel('Slowness', 'interpreter', 'latex', 'FontSize', 15);
ylabel('Depth index', 'interpreter', 'latex', 'FontSize', 15);

end

