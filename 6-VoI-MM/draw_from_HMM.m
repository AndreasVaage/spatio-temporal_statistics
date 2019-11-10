function x = draw_from_HMM(P,prio_x1,N)
%Draws a realisation of a binary first order Hiddem Markov Model
%   Input:
%   Markov transition matrix
%   P = [ p(x_{i+1} = 0 | x_i = 0), p(x_{i+1} = 0 | x_i = 1)
%         p(x_{i+1} = 1 | x_i = 0), p(x_{i+1} = 1 | x_i = 1)]
%
%   Prior x_1: 2x1 vector
%   prio_x1(x+1) = p(x_1 = x)         x = 0 or 1
%
%   Output
%   Realisation: x Nx1 vector
%   x(i) = 0 or 1

x = zeros(N,1);

if(rand<=prio_x1(1))
    x(1) = 0;
else
    x(1) = 1;
end

for k=2:N
    r = rand;
   if(x(k-1)==0)
       if(r<=P(1,1))
           x(k) = 0;
       else
           x(k) = 1;
       end
   else
       if(r<=P(2,2))
           x(k) = 1;
       else
           x(k) = 0;
       end
   end
end
end

