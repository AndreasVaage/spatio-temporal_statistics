clc
%close all
clear variables

%% Task 3

m = 200;            % Number of measurments
sigma2 = 1;         % Covariance parameter
eta = 10;           % Covariance parameter
tau2 = 0.05^2;      % Covariance parameter (Measurment Noise strddev)

alpha_true = 1;     % Regression parameter (Mean parameter)
theta_true = [sigma2; eta; tau2];

% Generate data sites

dataSites = rand(m,2); % data acquisition locations

% Generate Matern Covariance Matrix

Cov_true = CalcMaternCov(dataSites, theta_true);

% Generate measurements

C = chol(Cov_true)';
Y = C*randn(m,1);    % Y ~ N(0,Sigma)

H = zeros(m,1);
for i = 1:m
    H(i) = dataSites(i,1) + dataSites(i,2) - 1;
end

Y = Y + alpha_true * H;             % Y ~ N(alpha*H,Sigma)

%% Plotting

figure(1);
for i = 1:m
    plot(dataSites(i,1), dataSites(i,2), '.b'); hold on;
    text(dataSites(i,1), dataSites(i,2), [' ' num2str(Y(i))]);
end
hold off;
grid on; box on;
title(['\textbf{' num2str(m) ' Random Sample Locations}'], 'Interpreter', 'latex', 'FontSize', 15);
xlabel('East'); ylabel('North');
%%
figure(2);
imagesc(Cov_true);
title('\textbf{Matern Covariance Matrix}', 'Interpreter', 'latex', 'FontSize', 15);

%% Task 2

N_steps = 25;

theta_est = zeros(N_steps, 3);
%theta_est(1,:) = [3;5;0.001];  % Initial guess
theta_est(1,:) = [20;20;0.001];  % Initial guess
alpha_est = zeros(N_steps,1);
logLikelihood = zeros(N_steps,1);


for i = 1:(N_steps -1)
    
    Cov_est = CalcMaternCov(dataSites, theta_est(i,:));
    
    alpha_est(i) = (H'/Cov_est*H)\(H'/Cov_est*Y);
    Z = Y - H*alpha_est(i);
    
    
    logLikelihood(i) = -0.5*log(det(Cov_est)) - 0.5*Z'*(Cov_est\Z);
    
    
    [d_Cov_d_sig2, d_Cov_d_eta, d_Cov_d_tau2] = calc_d_Cov_d_theta(dataSites, theta_est(i,:));
    
    d_l_d_theta = zeros(3,1);
    
     % dl/dsigma2
    d_l_d_theta(1) = -0.5*trace(Cov_est\d_Cov_d_sig2) + 0.5*Z'*(Cov_est\d_Cov_d_sig2)*(Cov_est\Z);
    % dl/deta
    d_l_d_theta(2) = -0.5*trace(Cov_est\d_Cov_d_eta ) + 0.5*Z'*(Cov_est\d_Cov_d_eta )*(Cov_est\Z);
    % dl/dtau2
    d_l_d_theta(3) = -0.5*trace(Cov_est\d_Cov_d_tau2) + 0.5*Z'*(Cov_est\d_Cov_d_tau2)*(Cov_est\Z);
    
    expHessian_est = zeros(3,3);
    
    expHessian_est(1,1) = -trace((Cov_est\d_Cov_d_sig2)*(Cov_est\d_Cov_d_sig2));
    expHessian_est(1,2) = -trace((Cov_est\d_Cov_d_sig2)*(Cov_est\d_Cov_d_eta ));
    expHessian_est(1,3) = -trace((Cov_est\d_Cov_d_sig2)*(Cov_est\d_Cov_d_tau2));
    expHessian_est(2,2) = -trace((Cov_est\d_Cov_d_eta )*(Cov_est\d_Cov_d_eta ));
    expHessian_est(2,3) = -trace((Cov_est\d_Cov_d_eta )*(Cov_est\d_Cov_d_tau2));
    expHessian_est(3,3) = -trace((Cov_est\d_Cov_d_tau2)*(Cov_est\d_Cov_d_tau2));
    
    expHessian_est(2,1) = expHessian_est(1,2);
    expHessian_est(3,1) = expHessian_est(1,3);
    expHessian_est(3,2) = expHessian_est(2,3);


    theta_est(i+1,:) = theta_est(i,:) - (expHessian_est\d_l_d_theta)';
end

Cov_est = CalcMaternCov(dataSites, theta_est(i+1,:));
    
alpha_est(i+1) = (H'/Cov_est*H)\(H'/Cov_est*Y);
Z = Y - H*alpha_est(i+1);

logLikelihood(i+1) = -0.5*log(det(Cov_est)) - 0.5*Z'*(Cov_est\Z);


%% Task 3

gridX = 25;
gridY = 25;

n = gridX * gridY;

predictionSites = zeros(n, 2);
for i = 1:n
    predictionSites(i,1) = (2*mod(i-1, gridX)+1) / (2*gridX);
    predictionSites(i,2) = (2*floor((i-1)/gridX)+1) / (2*gridY);

end

figure(3);
for i = 1:n
    plot(predictionSites(i,1), predictionSites(i,2), '.b'); hold on;
    %text(siteCoords(i,1), siteCoords(i,2), [' ' num2str(Y(i))]);
end
hold off;
grid on; box on;
title(['\textbf{' num2str(n) ' Grid-Regular Prediction Sites}'], 'Interpreter', 'latex', 'FontSize', 15);
xlabel('East'); ylabel('North');

Cov_pred = CalcMaternCov(predictionSites, theta_est(N_steps,:));

%%

Y_pred = -ones(n,1);
Y_pred_var = -ones(n,1);

for i = 1:n
    [Y_pred(i), Y_pred_var(i)] = singlePointKrigingInterp(dataSites, predictionSites(i,:), theta_est(N_steps,:), Y);
end
%%
Y_pred_matrix = reshape(Y_pred,[25,25]);
Y_pred_var_matrix = reshape(Y_pred_var,[25,25]);
figure(4);
subplot(1,2,1);
heatmap(Y_pred_matrix);

subplot(1,2,2); hold off;
imagesc(Y_pred_var_matrix'); hold on;
for i = 1:m
    plot(25*dataSites(i,1), 25*dataSites(i,2), '.r'); hold on;
end
%heatmap(Y_pred_var_matrix); hold on;

%% Functions

function Cov = CalcMaternCov(siteCoords, theta)

m = length(siteCoords);
Cov = -1*ones(m,m);

for i = 1:m
    for j = 1:m
        dist_ij = norm(siteCoords(i,:) - siteCoords(j,:));
        Cov(i,j) = theta(1)*(1 + theta(2)*dist_ij) * exp(-theta(2)*dist_ij);
    end
end

end

function [val, val_var] = singlePointKrigingInterp(siteCoords, predCoords, theta, data)

cov = CalcMaternCov(siteCoords, theta);
cov_pred = CalcJointMaternCov(siteCoords, predCoords, theta);


val = (cov\cov_pred)'*data;
val_var = theta(3) - cov_pred'*(cov\cov_pred);


end

function Cov_pred = CalcJointMaternCov(siteCoords, predCoords, theta)

m = length(siteCoords);
n = size(predCoords);
n = n(1);

Cov_pred = -1*ones(m,n); % covariance of all the prediction sites to the data sites

for i = 1:m
    for j = 1:n
        dist_ij = norm(siteCoords(i,:) - predCoords(j,:));
        Cov_pred(i,j) = theta(1)*(1 + theta(2)*dist_ij) * exp(-theta(2)*dist_ij);
    end
end

end


function [d_Cov_d_sig2, d_Cov_d_eta, d_Cov_d_tau2] = calc_d_Cov_d_theta(siteCoords, theta_est)

m = length(siteCoords);

d_Cov_d_sig2 = -1*ones(m,m);
d_Cov_d_tau2 = -1*ones(m,m);
d_Cov_d_eta  = -1*ones(m,m);

for i = 1:m
    for j = 1:m
        dist_ij = norm(siteCoords(i,:) - siteCoords(j,:));
        
        d_Cov_d_sig2(i,j) = (1 + theta_est(2)*dist_ij) * exp(-theta_est(2)*dist_ij);
        
        d_Cov_d_eta(i,j) = - theta_est(2)*dist_ij^2 * exp(-theta_est(2)*dist_ij);
        
        d_Cov_d_tau2(i,j) = i == j;
    end
end

end

