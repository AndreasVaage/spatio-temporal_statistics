%% Second task.

samples = 1000;

rho_vec = -1:1/(samples/2):1;
tau_vec = [0.5,1];

voi1 = zeros(samples,1);
voi2 = zeros(samples,1);

% c)
for j=1:samples
    voi1(j) = voi(rho_vec(j),tau_vec(1));
    voi2(j) = voi(rho_vec(j),tau_vec(2));
end

figure
plot(rho_vec(2:end)',voi1);
hold on
plot(rho_vec(2:end)',voi2);

% d)
voi_tau_05 = zeros(samples,1);
voi_tau_1 = zeros(samples,1);
for i=1:samples
    voi_tau_05(i) = (1+abs(rho_vec(i)))/sqrt(2*pi*(1+tau_vec(1)^2));
    voi_tau_1(i) = (1+abs(rho_vec(i)))/sqrt(2*pi*(1+tau_vec(2)^2));
end

figure
plot(rho_vec(2:end)',voi_tau_05);
hold on
plot(rho_vec(2:end)',voi_tau_1);