function Sigma = generate_Sigma(theta, s, m)
sigma2 = theta(1);
eta = theta(2);
tau2 = theta(3);
C = zeros(m,m);
for i = 1:m
    for j = 1:m
        eucl_dist = norm(s(i)-s(j));
        C(i,j) = sigma2 * (1 + eta*eucl_dist) * exp(-eta * eucl_dist);
    end
end

Sigma = C + eye(m)*tau2;
end