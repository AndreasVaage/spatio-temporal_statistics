function marginal_likelihood_N = forward_reqursion(p,sigma,y,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
prio_x = [0.5,0.5];
P = [p 1-p; 1-p p]; % Markov transition matrix

like_yx = zeros(2,1);       % Likelihood
                            % like_yx(x_i) = p(y_i|x_i)
like_mean = [0;1];          % Mean of likelihood is x
like_std = sigma;  

prio_xY = zeros(2,N);       % Prior
                            % prio_xY(x,i) = p(x_i|y_1,...,y_{i-1})
p_yY = zeros(1,N);          % Normalising constant
                            % p_yY(i) = p(y_i|y_1,...,y_{i-1})
p_Y = zeros(1,N);           % Marginal likelihood 
                            % p_Y(i) = p(y_1,...y_i)

post_xY = zeros(2,N);       % Posteriori
                            % post_xY(x,i) = p(x_i|y_1,...,y_i)
post_mean = zeros(2,N);
post_std = zeros(2,N);

for i = 1:N

% Predicion
    if (i == 1)
        prio_xY(:,i)=prio_x;
    else
        prio_xY(:,i) = P*post_xY(:,i-1);
    end

% Filtering

    like_yx =  [normpdf(y(i),like_mean(1),like_std);
                normpdf(y(i),like_mean(2),like_std)];
    p_yY(i) = prio_xY(:,i)'*like_yx; % Matrix multiply to get sum for all x

    post_mean(:,i) = prio_xY(:,i).*like_mean./p_yY(i);
    post_std(:,i)  = prio_xY(:,i).*like_std./p_yY(i);
    post_xY(:,i) = [normpdf(0,post_mean(1,i),post_std(1,i));
                    normpdf(1,post_mean(2,i),post_std(2,i))];
    if (i == 1)
        p_Y(i) = p_yY(i);
    else
        p_Y(i) = p_yY(i)*p_Y(i-1);
    end
end
marginal_likelihood_N = p_Y(N);

end

