function marginal_log_likelihood_N = forward_reqursion(p,tau,y,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = [0;1];                  % x_i can either be 0 or 1
                            % For every probability the first row is true
                            % if x = 0, the second row if is true if x = 1
prio_x = [0.5;0.5];

P = [  p  1-p; 
     1-p    p];             % Markov transition matrix                            
                            
like_yx = zeros(2,N);       % Likelihood
                            % like_yx(x,i) = p(y_i|x_i)

like_yx(1,:) = normpdf(y, x(1), tau);
like_yx(2,:) = normpdf(y, x(2), tau);

prio_xY = zeros(2,N);       % Conditional Prior
                            % prio_xY(x,i) = p(x_i|y_1,...,y_{i-1})
p_yY = zeros(1,N);          % Normalising constant
                            % p_yY(i) = p(y_i|y_1,...,y_{i-1})
p_Y_log = zeros(1,N);       % Marginal likelihood 
                            % p_Y(i) = p(y_1,...y_i)

post_xY = zeros(2,N);       % Posteriori
                            % post_xY(x,i) = p(x_i|y_1,...,y_i)
%post_mean = zeros(2,N);
%post_std = zeros(2,N);

%%


%%
for i = 1:N

% Predicion
    if (i == 1)
        prio_xY(:,i) = prio_x;
    else
        prio_xY(:,i) = P*post_xY(:,i-1);
    end

% Filtering
        
    p_yY(i) = prio_xY(:,i)'*like_yx(:,i); % Matrix multiply to get sum for all x
    
    post_xY(:,i) = prio_xY(:,i) .* like_yx(:,i) ./ p_yY(i);
    
%    post_mean(:,i) = prio_xY(:,i) .* like_mean ./ p_yY(i);
%    post_std(:,i)  = prio_xY(:,i) .* like_std ./ p_yY(i);

% Marginalise
    if (i == 1)
        p_Y_log(i) = log(p_yY(i));
    else
        p_Y_log(i) = log(p_yY(i)) + p_Y_log(i-1);
    end
end
marginal_log_likelihood_N = p_Y_log(N);

end

