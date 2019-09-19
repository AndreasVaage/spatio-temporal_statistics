clear variables
clc
close all
N = 250;

s = 0.4*randn(N,1);
r = rand(N,1);

x = ones(N,1)*-1;
y = ones(N,1)*-1;

if r(1) < 0.5
    x(1) = 1;
else
    x(1) = 0;
end
%y(1) = x(1) + s(1);


for i = 2:N
    if r(i) < 0.9
        x(i) = x(i-1);
    else
        x(i) = 1 - x(i-1);
    end
    
    %y(i) = x(i) + s(i);
end

y = x + s;

%% Subplot plot
figure;
subplot(2,1,1);
plot(y); grid on;
ylabel('Data $y_i$', 'interpreter', 'latex', 'FontSize', 15);
xlabel('index $i$', 'interpreter', 'latex', 'FontSize', 15);
title('\textbf{Measurements $\mathbf y$}', 'interpreter', 'latex', 'FontSize', 15);

subplot(2,1,2);
stairs(x); grid on;
ylabel('States $x_i$', 'interpreter', 'latex', 'FontSize', 15);
xlabel('index $i$', 'interpreter', 'latex', 'FontSize', 15);
title('\textbf{Hidden States $\mathbf x$}', 'interpreter', 'latex', 'FontSize', 15);
axis([0 250 -0.5 1.5]);

%% Single graph plot
figure;
stairs(x); grid on; hold on;
plot(y);
ylabel('Data $y_i$, States $x_i$', 'interpreter', 'latex', 'FontSize', 15);
xlabel('index $i$', 'interpreter', 'latex', 'FontSize', 15);
title('\textbf{Hidden Markov Model Simulation}', 'interpreter', 'latex', 'FontSize', 15);
legend('States $\mathbf x$', 'Data $\mathbf y$', 'interpreter', 'latex');

%% B)
clc
close all
% Local
%taus = [0.3:0.005:0.5];
%ps = [0.85:0.001:0.999];
% Global
taus = [0.15:0.01:2];
ps = [0.2:0.01:0.999];

Z = zeros(length(ps),length(taus));

max = -inf;
r = 0;
for p = ps
    r = r + 1;
    c = 0;
    for tau = taus
        c = c + 1;
        
        Z(r,c) = forward_reqursion(p,tau,y,N);
        
        if (Z(r,c) >= max)
            max = Z(r,c);
            ml_p = p;
            ml_tau = tau;
        end
    end
end

disp("Max likelihood p = ");
disp(ml_p);
disp("Max likelihood tau = ");
disp(ml_tau);
%%
figure
[X,Y] = meshgrid(taus,ps);
mesh(X,Y,Z,'FaceAlpha',0.5,'FaceColor','interp')
xlabel("$\tau$",'interpreter', 'latex', 'FontSize', 15)
ylabel("$p$",'interpreter', 'latex', 'FontSize', 15)
zlabel("$\mathrm{log} \,\, p_{\theta}(\mathbf{y})$",'interpreter', 'latex', 'FontSize', 15)

%%
%C - Backward Recursion
clc
[~,prio_xY,post_xY] = forward_reqursion(0.9,0.4,y,N);
margp_xY = backward_recursion(0.9,prio_xY,post_xY,N);

disp("p(x_1 = 1 | Y) = "+ num2str(margp_xY(2,1)));

