function [ensamble, CovMat, L, mean] = genRealizations(B)

ensamble = -ones(100, B);

CovMat = -ones(100, 100);

sigma = 0.05;
mean = 0.5 - 0.001*(1:100)';
eta = 0.1;

for i = 1:100
    for j = 1:100
        
        CovMat(i, j) = sigma^2 * (1 + eta*abs(i - j))*exp(-eta*abs(i - j));
        
    end
end

L = chol(CovMat)';

for i = 1:B
    ensamble(:, i) = mean + L*randn(100, 1);
end


end