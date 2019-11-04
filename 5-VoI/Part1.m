close all
clear

directory = 'plots/part1/';

PoV = @(mu,sigma) mu .* normcdf(mu./sigma) + sigma .* normpdf(mu./sigma);

% i) Range of mu
mu_plt = -2:0.1:2;
% ii) Range of sigma
sigma_plt = 0.1:0.1:2;
[MU_plt,SIGMA_plt] = meshgrid(mu_plt, sigma_plt);
[r,c] = size(SIGMA_plt);
PoV_plt = zeros(r,c);
for i=1:1:r
for j=1:1:c
PoV_plt(i,j) = PoV(MU_plt(i,j),SIGMA_plt(i,j))-max(0,MU_plt(i,j));
end
end

fig1 = figure;
surf(MU_plt,SIGMA_plt,PoV_plt)
grid('on')
xlabel('$\mu$ ','Interpreter','latex', 'FontSize', 11)
ylabel('$\sigma$ ','Interpreter','latex', 'FontSize', 11)
zlabel('$PoV$','Interpreter','latex', 'FontSize', 11)
legend('Interpreter','latex', 'FontSize', 11)
% title_pic = 'PoV';
% save_plots(title_pic,directory,fig1);

%%
sigma_plt2 = [0.1,1,2];
mu_plt2 = -2:0.01:2;
fig2 = figure;
subplot(1,2,1)
hold on
for j = 1:length(sigma_plt2)
plot(mu_plt2,PoV(mu_plt2,sigma_plt2(j))-max(0,mu_plt2),'color',[0 0 0]+(j-1)*0.3,'DisplayName',strcat('$\sigma=',num2str(sigma_plt2(j)),'$'));
legend('Interpreter','latex', 'FontSize', 11)
xlabel('$\mu$','Interpreter','latex', 'FontSize', 11)
ylabel('VoI')
end
%
subplot(1,2,2)
hold on
sigma_plt3 = 0.01:0.01:2;
mu_plt3 = [-1,0,2];
hold on
for j = 1:length(mu_plt3)
plot(sigma_plt3,PoV(mu_plt3(j),sigma_plt3)-max(0,mu_plt3(j)),'color',[0 0 0]+(j-1)*0.3,'DisplayName',strcat('$\mu=',num2str(mu_plt3(j)),'$'));
legend('Interpreter','latex', 'FontSize', 11)
xlabel('$\sigma$','Interpreter','latex', 'FontSize', 11)
ylabel('VoI')
end
% title_pic = 'figg2';
% save_plots(title_pic,directory,fig2);
% myeplot