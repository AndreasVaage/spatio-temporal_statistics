function [ensamble, CovMat, L, mean_] = genRealizations(B)

layers = 100;

ensamble = -ones(layers, B);

CovMat = -ones(layers, layers);

sigma = 0.05;
mean_ = 0.5 - 0.001*(1:layers)';
eta = 0.1;

for i = 1:layers
    for j = 1:layers
        
        CovMat(i, j) = sigma^2 * (1 + eta*abs(i - j))*exp(-eta*abs(i - j));
        
    end
end

L = chol(CovMat)';

for i = 1:B
    ensamble(:, i) = mean_ + L*randn(layers, 1);
end


end