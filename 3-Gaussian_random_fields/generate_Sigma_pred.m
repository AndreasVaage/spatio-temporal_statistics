function C_0 = generate_Sigma_pred(theta, s_pred, s_meas, m)
sigma2 = theta(1);
eta = theta(2);
tau2 = theta(3);
C_0 = zeros(1,m);
for j = 1:m
    eucl_dist = norm( s_pred - s_meas(j,:) ); %sqrt((s_pred(1)-s_meas(j,1))^2 + (s_pred(2)-s_meas(j,2))^2)
    C_0(1,j) = sigma2 * (1 + eta*eucl_dist) * exp(-eta * eucl_dist);
end

end