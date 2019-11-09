function margp_xY = backward_recursion(p,prio_xY,post_xY,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


P = [  p  1-p;
       0    1]';        % Markov transition matrix                            
                            
%prio_xY = zeros(2,N); % Conditional Prior
                       % prio_xY(x,i) = p(x_i | y_1,...,y_{i-1})
                            
margp_xY = zeros(2,N); % Marginal probabilities
                       % margp_xY(x,i) = p(x_i | y_1,...,y_n)
p_xYx = zeros(2,2);    % p_xYx(x_i,x_{i+1}) = p(x_i|y_1,...,y_{i+1},x_{i+1}
                            
%post_xY = zeros(2,N); % Posteriori
                       % post_xY(x,i) = p(x_i|y_1,...,y_i)

%%

margp_xY(:,N) = post_xY(:,N);

for i = N-1:-1:1
    
    %p_xYx = P' .* post_xY(:,i) ./ prio_xY(:,i+1);
    
    p_xYx(1,1) = P(1,1) * post_xY(1,i) / prio_xY(1,i+1);
    p_xYx(1,2) = P(2,1) * post_xY(1,i) / prio_xY(2,i+1);
    
    p_xYx(2,1) = P(1,2) * post_xY(2,i) / prio_xY(1,i+1);
    p_xYx(2,2) = P(2,2) * post_xY(2,i) / prio_xY(2,i+1);
    
    margp_xY(:,i) = p_xYx(:,:)*margp_xY(:,i+1);

end

end

