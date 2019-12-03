function [dLogLikelihood, Hessian] = analythic_derivation(Sigma, Z,theta, s, m)
sigma2 = theta(1);
eta = theta(2);
%tau2 = theta(3);

dC_eta = zeros(m,m);
dC_sigma2 = zeros(m,m);
for i = 1:m
    for j = 1:m
        eucl_dist = norm(s(i)-s(j));
        dC_sigma2(i,j) = (1 + eta*eucl_dist)*exp(-eta * eucl_dist);
        dC_eta(i,j) = sigma2 *( -eta*(eucl_dist^2)*exp(-eta * eucl_dist));
    end
end
dSigma_eta = dC_eta;
dSigma_sigma2 = dC_sigma2;
dSigma_tau2 = eye(m);

%theta = zeros(3,1); %[sigma^2 eta tau^2]'

%Q = inv(Sigma);

Hessian = zeros(3); % Expected value of the Hessian of the loglikelyhood

QdSigma_sigma2 = Sigma\dSigma_sigma2;
QdSigma_eta = Sigma\dSigma_eta;
QdSigma_tau2 = Sigma\dSigma_tau2;
Hessian(1,1) = -(1/2)*trace(QdSigma_sigma2*QdSigma_sigma2);
Hessian(1,2) = -(1/2)*trace(QdSigma_sigma2*QdSigma_eta);
Hessian(1,3) = -(1/2)*trace(QdSigma_sigma2*QdSigma_tau2);
Hessian(2,1) = Hessian(1,2);
Hessian(2,2) = -(1/2)*trace(QdSigma_eta*QdSigma_eta);
Hessian(2,3) = -(1/2)*trace(QdSigma_eta*QdSigma_tau2);
Hessian(3,1) = Hessian(1,3);
Hessian(3,2) = Hessian(2,3);
Hessian(3,3) = -(1/2)*trace(QdSigma_tau2*QdSigma_tau2);

dLogLikelihood = zeros(3,1);
%ZQ = Z'/Sigma;
dLogLikelihood(1) = -(1/2)*trace(QdSigma_sigma2) + (1/2)*Z'*(QdSigma_sigma2/Sigma)*Z;
dLogLikelihood(2) = -(1/2)*trace(QdSigma_eta) + (1/2)*Z'*(QdSigma_eta/Sigma)*Z;
dLogLikelihood(3) = -(1/2)*trace(QdSigma_tau2) + (1/2)*Z'*(QdSigma_tau2/Sigma)*Z;
end