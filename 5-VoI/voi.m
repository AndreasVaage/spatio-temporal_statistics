function v = voi(rho,tau)
%POV Summary of this function goes here
%   Detailed explanation goes here

sigma = [1 rho; rho 1];

cov = sigma*((sigma + eye(2)*tau)\sigma); %inv(sigma + eye(2)*tau)*sigma

v = sum(sqrt(diag(cov)))/sqrt(2*pi);

end

