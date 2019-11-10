function [marginal_log_likelihood_N, prio_xY, post_xY] = forward_recursion(P,prio_x1,tau,y,N)
%Calculates the posterior probability using forward recursion
%   Input:
%   Markov transition matrix
%   P = [ p(x_{i+1} = 0 | x_i = 0), p(x_{i+1} = 0 | x_i = 1)
%         p(x_{i+1} = 1 | x_i = 0), p(x_{i+1} = 1 | x_i = 1)]
%
%   Prior x_1: prio_x1 2x1 matrix
%   prio_x1(x+1) = p(x_1 = x)                   x = 0 or 1
%
%   Measurments: y 1xN vector
%   Unbiased conditional independetn meausurments with
%   Measurment standard deviation: tau 1xN vector or scalar
%   If tau is a scalar then it is assumend all measurments have the same
%   standard deviation
%   y(i) ~ N(x(i),tau(i)^2), i = 1:N
%
%   Output:
%   Marginal log likelihood
%   log(p(y_1, ..., y_n))   Probability of getting those exact measurments
%
%   Conditional Prior: Prio_xY 2xN matrix
%   prio_xY(x+1,i) = p(x_i = x|y_1,...,y_{i-1}) x = 0 or 1
%
%   Posteriori: Post_xY 2xN matrix
%   post_xY(x+1,i) = p(x_i = x|y_1,...,y_i)     x = 0 or 1


x = [0;1];                  % x_i can either be 0 or 1
                            % For every probability the first row is true
                            % if x = 0, the second row if is true if x = 1
                            
like_yx = zeros(2,N);       % Likelihood
                            % like_yx(x,i) = p(y_i|x_i)

assert( isequal(size(tau),size(y)) || isequal(size(tau),[1,1]),...
    "tau has to be same size as y or a scalar")
 
like_yx(1,:) = normpdf(y, x(1,:), tau);
like_yx(2,:) = normpdf(y, x(2,:), tau);
                         

prio_xY = zeros(2,N);       % Conditional Prior
                            % prio_xY(x,i) = p(x_i|y_1,...,y_{i-1})
p_yY = zeros(1,1);          % Normalising constant
                            % p_yY = p(y_i|y_1,...,y_{i-1})
p_Y_log = zeros(1,N);       % Marginal likelihood 
                            % p_Y(i) = p(y_1,...y_i)

post_xY = zeros(2,N);       % Posteriori
                            % post_xY(x,i) = p(x_i|y_1,...,y_i)
%post_mean = zeros(2,N);
%post_std = zeros(2,N);

%%
for i = 1:N

% Predicion
    if (i == 1)
        prio_xY(:,i) = prio_x1;
    else
        prio_xY(:,i) = P*post_xY(:,i-1);
    end

% Filtering
    
    
%%
    p_yY = prio_xY(:,i)'*like_yx(:,i);
    
    
    post_xY(:,i) = prio_xY(:,i) .* like_yx(:,i) ./ p_yY;


    if (i == 1)
        p_Y_log(i) = log(p_yY);
    else
        p_Y_log(i) = log(p_yY) + p_Y_log(i-1);
    end
end

marginal_log_likelihood_N = p_Y_log(N);

% size(prio_xY)
% size(post_xY)

end

