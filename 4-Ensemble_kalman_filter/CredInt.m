function [yMean, yCIpercen] = CredInt(y,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by: J. Mendoza
% Date: January, 2019
% E: jorge.m.espinosa@ntnu.no
% Insitution: Norwegian University of Science and Technology (NTNU)
% Place: Trondheim, Norway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Confidence interval of the experiments at each value of x
% x : Independent Variable [1xM] array
% INPUT
% y : Dependent Variable ‘Experiments’ Data. [NXM]
% percen: 0.95 by default.
% OUTPUT
% yMean: % Mean of y for all experiments at each value of x
%
switch nargin
    case 1
       percen = 0.95;
    case 2
       percen = varargin{1};
end
p = [0.5*(1-percen) 0.5*(1+percen)]; 

N = size(y,1);                                                             % Number of experiments in data set
yMean = mean(y,1);                                                         % Mean of y for all experiments at each value of x
ySEM = std(y)/sqrt(N);                                                     % Compute standard error of the mean of all experiments at each value of x
CIp = tinv(p, N-1);                                                        % Calculate p probability intervals of t-Distribution
yCIpercen = bsxfun(@times, ySEM, CIp(:));                                  % Calculate p confidence intervals of all experiments at each value of ‘x’