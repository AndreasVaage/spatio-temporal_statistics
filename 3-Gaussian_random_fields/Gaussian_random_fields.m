%% Gaussian processes
% Here we will study Gaussian processes or random fields. We simulate spatial 
% data, perform parameter estimation and prediction (Kriging).
% 
% 
%% Task 1
% 
clc
close all
clear variables

m = 200;       % Number of measurments
sigma2 = 1;    % Covariance parameter
eta = 10;      % Covariance parameter
tau2 = 0.05^2; % Covariance parameter (Measurment Noise std)
alpha_true = 1;     % Regression parameter (Mean parameter)
theta_true = [sigma2; eta; tau2];

s = rand(200,2);
figure
plot(s(:,1),s(:,2),'.');

Sigma = generate_Sigma(theta_true, s, m);

%Y = chol(C)*rand(m,1) + tau*randn(m,1);

Y = chol(Sigma)'*randn(m,1);    % Y ~ N(0,Sigma)

H = zeros(m,1);
for i = 1:m
    H(i) = (s(i,1) - 0.5) + (s(i,2) - 0.5);
end

Y = Y + alpha_true * H;             % Y ~ N(alpha*H,Sigma)
%% Task 2
% 
% 
% Log likelihood of multivariate gaussian:
% 
% $$l\left(\mathbf{Y}\;\left|\;\alpha ,\theta \;\right.\right)=-\frac{1}{2}\log 
% \;\left|\Sigma \right|-\frac{1}{2}\;{\mathbf{Z}}^{\prime } {\;\Sigma }^{-1} 
% \;\mathbf{Z}$$
% 
% $$\mathbf{Z}=\mathbf{Y}-\mathbf{H}\alpha$$
% 
% $$\theta =\left\lbrack \begin{array}{ccc}\sigma^2  & \tau^2  & \eta \end{array}\right\rbrack$$
% 
% 
% 
% Derivate loglikelihood with respect to $\theta$:
% 
% 
% 
% 
% 
% Derivate loglikelihood with respect to $\alpha$:
% 
% $$\frac{d}{\mathrm{d}\alpha }l=-{\mathbf{H}}^{\prime } \;\mathbf{Q}\;\mathbf{H}\;\alpha 
% +{\mathbf{H}}^{\prime } \;\mathbf{Q}\;\mathbf{Y}$$
% 
% Set to zero:
% 
% $$\begin{array}{l}\alpha^ˆ ={{\left({\mathbf{H}}^{\prime } \;\mathbf{Q}\;\mathbf{H}\;\right)}^{-1} 
% \;{\mathbf{H}}^{\prime } \;\mathbf{Q}\;\mathbf{Y}\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;}^{\prime 
% } \\\alpha^ˆ ={\mathbf{H}}^{-1} \;\mathbf{Y}\end{array}$$
% 
% NB: This estimated $\alpha \;$is a function of $\theta$ thus we update them 
% iteratively.
% 
% $\beta =\alpha$

N_steps = 10;
theta = zeros(3,N_steps);
theta(:,1) = [3;5;0.0001];  % Initial guess
alpha = zeros(N_steps,1);
LogLikelihood = zeros(N_steps,1);
for i = 1:(N_steps -1)
    Sigma = generate_Sigma(theta(:,i), s, m);
    %Q = Sigma\eye(m);
    %alpha(i) = H\Y 
    alpha(i) = (H'/Sigma*H)\(H'/Sigma*Y);
    %alpha(i) = (H'/Sigma*H)\(H'/Sigma*Y)
    Z = Y - H*alpha(i);
    LogLikelihood(i) = -0.5*log(norm(Sigma)) - 0.5*Z'/Sigma*Z;
    
    [dLogLikelihood, Hessian] = analythic_derivation(Sigma, Z,theta(:,i), s, m);
    
    theta(:,i+1) = theta(:,i) - Hessian\dLogLikelihood;
end
    Sigma = generate_Sigma(theta(:,i+1), s, m);
    alpha(i+1) = (H'/Sigma*H)\(H'/Sigma*Y);
    Z = Y - H*alpha(i+1);
    LogLikelihood(i+1) = -0.5*log(norm(Sigma)) - 0.5*Z'/Sigma*Z;
    
figure
subplot(5,1,1); hold on;
yline(alpha_true,'-.');
plot(alpha)
legend("$\hat{\alpha}$","$\alpha$",'interpreter', 'latex', 'FontSize', 15)
subplot(5,1,2); hold on;
yline(theta_true(1),'-.');
plot(theta(1,:))
legend("$\hat{\sigma}^2$","$\sigma^2$",'interpreter', 'latex', 'FontSize', 15)
subplot(5,1,3); hold on;
yline(theta_true(2),'-.');
plot(theta(2,:))
legend("$\hat{\eta}$","$\eta$",'interpreter', 'latex', 'FontSize', 15)
subplot(5,1,4); hold on;
yline(theta_true(3),'-.');
plot(theta(3,:))
legend("$\hat{\tau}^2$","$\tau^2$",'interpreter', 'latex', 'FontSize', 15)
subplot(5,1,5); hold on;
plot((LogLikelihood))
legend("$l(\mathbf{Y} | \hat\theta,\hat\alpha)$",'interpreter', 'latex', 'FontSize', 15)
%% Task 3
% 
% 
% 
% 
% 
% 
%