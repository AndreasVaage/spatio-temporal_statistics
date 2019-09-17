
n = 250;

s = 0.4*randn(n,1);
r = rand(n,1);

x = ones(n,1)*-1;
y = ones(n,1)*-1;

if r(1) < 0.5
    x(1) = 1;
else
    x(1) = 0;
end
%y(1) = x(1) + s(1);


for i = 2:n
    if r(i) < 0.9
        x(i) = x(i-1);
    else
        x(i) = 1 - x(i-1);
    end
    
    %y(i) = x(i) + s(i);
end

y = x + s;

%%
figure; plot(y); grid on;
hold on; stem(x);
title('\textbf{Binary Hidden Markov Chain Simulation}', 'interpreter', 'latex', 'FontSize', 18);
xlabel('index $i$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('Data $\mathbf y$', 'interpreter', 'latex', 'FontSize', 15);
legend('Data y');

%% B)
%Init
sigma = 0.4;
py1= normpdf(y(1),0,sigma)*0.5 + normpdf(y(1),1,sigma)*0.5;

p = zeros(250,1);
p(1) = normpdf(y(1),0,sigma)*0.5 + normpdf(y(1),1,sigma)*0.5;

for i=2:250
    p(i) = (normpdf(y(i-1),0,sigma)*0.5 + normpdf(y(i-1),1,sigma)*0.5)*p(i-1);

p_filter_1 = 0.5*[x1,sigma]/py1;

filter_sigma = 0.5*sigma/py1;
filter_mean = 0.5/py1;




