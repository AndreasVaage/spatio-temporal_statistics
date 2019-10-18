PoV = @(mu,sigma) mu .* normcdf(mu/sigma) + sigma .* normpdf(mu/sigma);

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

figure(1)
surf(MU_plt,SIGMA_plt,PoV_plt)
grid('on')
xlabel('$\mu$ ')
ylabel('$\sigma$ ')
zlabel('$PoV$')
myeplot
%%
sigma_plt2 = [0.1,1,2];
mu_plt2 = -2:0.01:2;
figure
hold on
for j = 1:length(sigma_plt2)
plot(mu_plt2,PoV(mu_plt2,sigma_plt2(j))-max(0,mu_plt2),'color',[0 0 0]+(j-1)*0.3,'DisplayName',strcat('$\sigma=',num2str(sigma_plt2(j)),'$'));
legend('interpreter','latex')
xlabel('$\mu$')
ylabel('VoI')
end