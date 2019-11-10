function margp_xY = backward_recursion(P,prio_xY,post_xY,N)
%Based on the output from forward recursion, calculates the total
%conditional propability p(x_i = x | y_1,...,y_n)
%   Input:
%   Markov transition matrix
%   P = [ p(x_{i+1} = 0 | x_i = 0), p(x_{i+1} = 0 | x_i = 1);
%         p(x_{i+1} = 1 | x_i = 0), p(x_{i+1} = 1 | x_i = 1)]
%
%   Conditional Prior prio_xY 2xN matrix
%   prio_xY(x+1,i) = p(x_i = x|y_1,...,y_{i-1}) x = 0 or 1
%
%   Posteriori Post_xY 2xN matrix
%   post_xY(x+1,i) = p(x_i = x | y_1,...,y_i)   x = 0 or 1
%
%   Output:
%   Marginal probabilities margp_xY 2xN matrix
%   margp_xY(x+1,i) = p(x_i = x | y_1,...,y_n)  x = 0 or 1
                      
                            
margp_xY = zeros(2,N); % Marginal probabilities
                       % margp_xY(x,i) = p(x_i | y_1,...,y_n)
                       
p_xYx = zeros(2,2);    % p_xYx(k+1,x_{i+1}) = p(x_i = k | y_1,...,y_{i+1},x_{i+1}


%%

margp_xY(:,N) = post_xY(:,N);

for i = N-1:-1:1
    
    %p_xYx = P' .* post_xY(:,i) ./ prio_xY(:,i+1);
    
    % p_xYx(1,2) = p(x_i = 0 | y_1,...,y_{i+1},x_{i+1})
    %            = p(x_{i+1} = 1|x_i = 0) * p(x_i = k | y_1,...,y_i)/ p(x_i = 0|y_1,...,y_{i-1})
    
    p_xYx(1,1) = P(1,1) * post_xY(1,i) / prio_xY(1,i+1);
    p_xYx(1,2) = P(2,1) * post_xY(1,i) / prio_xY(2,i+1);
    
    p_xYx(2,1) = P(1,2) * post_xY(2,i) / prio_xY(1,i+1);
    p_xYx(2,2) = P(2,2) * post_xY(2,i) / prio_xY(2,i+1);
    
    margp_xY(:,i) = p_xYx(:,:)*margp_xY(:,i+1);

end

end

