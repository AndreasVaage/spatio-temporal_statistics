%% Second task.

close all
clear

directory = 'plots/part2/';

samples = 1000;

rho_vec = -1:1/(samples/2):1;
tau_vec = [0.5,1];

voi1 = zeros(samples,1);
voi2 = zeros(samples,1);

% c)
for j=1:samples
    voi1(j) = voi(rho_vec(j),tau_vec(1));
    voi2(j) = voi(rho_vec(j),tau_vec(2));
end

fig1 = figure;
plot(rho_vec(2:end)',voi1);
hold on
plot(rho_vec(2:end)',voi2);
xlabel('$\rho$ ','Interpreter','latex', 'FontSize', 11)
ylabel('$VoI$ ','Interpreter','latex', 'FontSize', 11)
legend('$\tau=0.5$','$\tau=1$','Interpreter','latex', 'FontSize', 13);
% title_pic = 'VoI';
% save_plots(title_pic,directory,fig1);

% d)
voi_tau_05 = zeros(samples,1);
voi_tau_1 = zeros(samples,1);
for i=1:samples
    voi_tau_05(i) = (1+abs(rho_vec(i)))/sqrt(2*pi*(1+tau_vec(1)^2));
    voi_tau_1(i) = (1+abs(rho_vec(i)))/sqrt(2*pi*(1+tau_vec(2)^2));
end

fig2 = figure;
plot(rho_vec(2:end)',voi_tau_05);
hold on
plot(rho_vec(2:end)',voi_tau_1);
xlabel('$\rho$ ','Interpreter','latex', 'FontSize', 11)
ylabel('$VoI$ ','Interpreter','latex', 'FontSize', 11)
legend('$\tau=0.5$','$\tau=1$','Interpreter','latex', 'FontSize', 13);
% title_pic = 'VoI2';
% save_plots(title_pic,directory,fig2);