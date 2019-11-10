clc
%clear variables
close all

GS = 25;
N = GS^2;

s = -ones(N,2);
for r = 1:GS
    i = (r-1)*GS + 1;
    s(i:i+GS-1,1) = r*ones(GS,1);
    s(i:i+GS-1,2) = 1:GS;
end

Sigma = -ones(N,N);
for i = 1:N
    for j = 1:N
        Sigma(i,j) = exp(-0.15*norm(s(i,:)-s(j,:)));
    end
end

T = diag(5^2*ones(1,GS));

%% Part a) - Variance minimisation

Sigma_j = Sigma;
selected_rows = zeros(1,GS);
for t = 1:10
    j = varMinimizingRow(Sigma_j,T,GS);
    %plot2dData(diag(Sigma_j),GS);
    Sigma_j = updateCovariance(j,Sigma_j,T,GS);
    selected_rows(j) = t;
end
disp("Rows numbered by order of selection:");
disp(selected_rows);
%image(Sigma,'CDataMapping','scaled')


%% Part b) - VOI

%x = chol(Sigma)'*randn(N,1);
plot2dData(x,GS);

mu = zeros(N,1);
P = 0.5;

VOIs = zeros(GS,1);
for j = 1:GS
    VOIs(j) = calculateVOI(j,Sigma,mu,T,GS);
end
figure
plot(VOIs,'*')
xlabel("Data Sampeling Line, j")
ylabel("VOI_j")
%%
n_selected_rows = zeros(100,GS);
for game_number = 1:100
    Sigma_j = Sigma;
    mu_j = mu;
    selected_rows = zeros(1,GS);
    [VOI,j] = VOIMaximizinRow(Sigma_j,mu_j,T,GS);
    while(VOI > P)
        selected_rows(j) = selected_rows(j) + 1;
        y_j = sampleData(j,x,T,GS);
        mu_j = updateMean(j,Sigma_j,mu_j,y_j,T,GS);
        Sigma_j = updateCovariance(j,Sigma_j,T,GS);
        [VOI,j] = VOIMaximizinRow(Sigma_j,mu_j,T,GS);
    end
    n_selected_rows(game_number,:) = selected_rows;
end
figure
histogram(sum(n_selected_rows,2));
xlabel("Number of selected lines")
%%
A = [sum(n_selected_rows,2),n_selected_rows];
sorted_selected_rows = sortrows(A);
figure
image(sorted_selected_rows(:,2:end),'CDataMapping','scaled')
%set(gca,'YDir','normal')
colorbar
ylabel ("Trial");
xlabel("Data Sampeling Line, j")


%%
function plot2dData(x,GS)
    X = reshape(x,GS,GS)';
    figure
    image(X','CDataMapping','scaled')
    set(gca,'YDir','normal')
    colorbar
    ylabel ("North");
    xlabel("West");
end

function F = generateF(j,GS)
assert(j<=GS && j>0,"Index out of range");
F = zeros(GS,GS^2);
F(:,(j-1)*GS + 1:(j-1)*GS + GS) = eye(GS);
end

function y = sampleData(j,x,T,GS)
% Sample data alog row j, from x with Covarianse T of size GS
F = generateF(j,GS);
y = F*x + chol(T)'*randn(GS,1);
end

function Sigma_j = updateCovariance(j,Sigma,T,GS)
    F = generateF(j,GS);
    Sigma_j = Sigma - (Sigma*F')/(F*Sigma*F' + T)*(F*Sigma);
end

function mu_j = updateMean(j,Sigma,mu,y,T,GS)
    F = generateF(j,GS);
    mu_j = mu + (Sigma*F')/(F*Sigma*F' + T)*(y - F*mu);
end

function j = varMinimizingRow(Sigma,T,GS)
% Given current Sigma and T returns the variance minimizing row for
% sampeling data in a grid with gridsize GS
mean_variance = zeros(GS,1);
for j = 1:GS
    Sigma_j = updateCovariance(j,Sigma,T,GS);
    mean_variance(j) = mean(diag(Sigma_j));
end
[~,j] = min(mean_variance);
end

function VOI = calculateVOI(j,Sigma,mu,T,GS)
    assert(j<=GS && j>0,"Index out of range");
    assert(size(Sigma,2) == size(mu,1) && size(mu,2) == 1,"Wrong matrix dimention");
    mu_w = sum(mu);
    F_j = generateF(j,GS);
    R_j = (Sigma*F_j')/(F_j*Sigma*F_j' + T)*(F_j*Sigma);
    r_wj = sqrt(sum(sum(R_j)));
    PV = max(0,mu_w);
    PoV = mu_w * normcdf(mu_w/r_wj) + r_wj * normpdf(mu_w/r_wj);
    VOI = PoV - PV;
end

function [VOI,j] = VOIMaximizinRow(Sigma,mu,T,GS)
% Given current Sigma, mu and T returns the VOI maximizing row for
% sampeling data in a grid with gridsize GS
VOIs = zeros(GS,1);
for j = 1:GS
    VOIs(j) = calculateVOI(j,Sigma,mu,T,GS);
end
[VOI,j] = max(VOIs);
end

