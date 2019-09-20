function marginal_log_likelihood_N = forward_recursion(p,sigma,y,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


P0 = p;
P1 = 1-p;

p0 = 0.5;
p1 = 1-p0;


l0 = zeros(1,N);
l1 = zeros(1,N);

post0 = zeros(1,N);
post1 = zeros(1,N);

pred0 = zeros(1,N);
pred1 = zeros(1,N);

log_likelihood = 0;

for i = 1:N
    l0(i) = p0*normpdf(y(i),0,sigma);
    l1(i) = p1*normpdf(y(i),1,sigma);
    
    post0(i) = l0(i)/(l0(i) + l1(i));
    post1(i) = l1(i)/(l0(i) + l1(i));
    
    log_likelihood = log_likelihood + log(l0(i) + l1(i));
    
    pred0(i) = P0*post0(i) + P1*post1(i);
    pred1(i) = 1- pred0(i);
    
    p0 = pred0(i);
    p1 = pred1(i);
end

marginal_log_likelihood_N = log_likelihood;

end

