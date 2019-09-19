
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


for it_p = 2:N
    if r(it_p) < 0.9
        x(it_p) = x(it_p-1);
    else
        x(it_p) = 1 - x(it_p-1);
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
%Init
% sigma = 0.4;
% py1= normpdf(y(1),0,sigma)*0.5 + normpdf(y(1),1,sigma)*0.5;
% 
% p = zeros(250,1);
% p(1) = normpdf(y(1),0,sigma)*0.5 + normpdf(y(1),1,sigma)*0.5;
% 
% for i=2:250
%     p(i) = (normpdf(y(i-1),0,sigma)*0.5 + normpdf(y(i-1),1,sigma)*0.5)*p(i-1);
% 
% p_filter_1 = 0.5*[x1,sigma]/py1;
% 
% filter_sigma = 0.5*sigma/py1;
% filter_mean = 0.5/py1;
%%
clc
close all
% Local
sigmas = [0.3:0.005:0.6];
ps = [0.85:0.005:0.93];
% Global
%sigmas = [0:0.01:1];
%ps = [0.01:0.01:0.99];

Z = zeros(length(ps),length(sigmas));

it_p = 1;
max = 0;

ml_p = 0;
ml_sigma = 0;

for p = ps
    
    it_sigma = 1;
    
    for sigma = sigmas
        
        Z(it_p,it_sigma) = forward_reqursion(p,sigma,y,N);
        
        if (Z(it_p,it_sigma) >= max)
            max = Z(it_p,it_sigma);
            ml_p = p;
            ml_sigma = sigma;
        end
        it_sigma = it_sigma+1;
    end
    it_p = it_p+1;
end

sigma_ml =  ml_sigma
p_ml = ml_p

figure
[X,Y] = meshgrid(sigmas,ps);
%C = log(Z);
mesh(X,Y,Z,'FaceAlpha',0.5,'FaceColor','interp')
%set(gca,'ZScale','log')
%set(gca,'ColorScale','log')
xlabel("$\tau$",'interpreter', 'latex', 'FontSize', 15)
ylabel("$p$",'interpreter', 'latex', 'FontSize', 15)
zlabel("$p_{\theta}(\mathbf{y})$",'interpreter', 'latex', 'FontSize', 15)



