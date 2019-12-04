close all
clear variables
clc

directory_a = 'plots/a/';
directory_b = 'plots/b/';
directory_c = 'plots/c/';
directory_d = 'plots/d/';
directory_e = 'plots/e/';
directory_f = 'plots/f/';

%% task a
N = 50;
p = 0.95;       % p = p(x_{i+1} = 0 | x_i = 0)
q = 1;          % q = p(x_{i+1} = 1 | x_i = 1)
P = [p, 1-q;    % P = [ p(x_{i+1} = 0 | x_i = 0), p(x_{i+1} = 0 | x_i = 1)
    1-p, q];    %       p(x_{i+1} = 1 | x_i = 0), p(x_{i+1} = 1 | x_i = 1)]
% Note: P := Lecture_P'
marg_p = zeros(N,1);
marg_p(1) = 0.01;
p_x = zeros(2,N);   % p_x(x+1,k) = p(x_k = x)        x = 0 or 1
p_x(:,1) = [0.99;0.01];
for k=2:N
    marg_p(k)= marg_p(k-1) + 0.05*(1-marg_p(k-1));
    p_x(:,k) = P * p_x(:,k-1);
end
figure(1);hold off
plot(p_x(2,:));hold on
plot(p_x(1,:))
xlabel('railroad section $i$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$p(x_i)$', 'interpreter', 'latex', 'FontSize', 16);
legend('$p(x_i = 1)$','$p(x_i = 0)$', 'Location', 'best', 'interpreter', 'latex', 'FontSize', 14);
%% task b

[PV,a] = max([-100000, -5000*sum(p_x(2,:))]);
PV
a = a-1

%% task c
tau = 0.3;
x = draw_from_HMM(P,p_x(:,1),N);

y = normrnd(x,tau,N,1); % same as mu + sigma*randn(data,1)
y(20) = 0.2;
y(30) = 0.7;
v_tau = 100*tau*ones(size(y));
v_tau(20) = tau;
v_tau(30) = tau;

[~,prio_xY,post_xY] = forward_recursion(P,p_x(:,1),v_tau,y,N);
margp_xY = backward_recursion(P,prio_xY,post_xY,N);

figure(2);hold off
plot(1:N,margp_xY(2,:));
hold on
plot(1:N,post_xY(2,:),'--');
plot([20,30],[0.2,0.7],'*');
grid
xlabel('railroad section $i$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('probability', 'interpreter', 'latex', 'FontSize', 16);
legend('$p(x_i = 1|\mathbf{y}_D)$','$\mathbf{y}_D$', 'Location', 'best', 'interpreter', 'latex', 'FontSize', 14);
%save_plots(title_pic,directory_c,fig3);

%% task d
value_clean_all = -100000;
PV = -100000;

% Monte carlo
PoV = -ones(N,1);
mc_samples = 10000;
for k=1:N
    mean_PoV = 0; 
    for j=1:mc_samples
        x = draw_from_HMM(P,p_x(:,1),N);
        y = normrnd(x,tau,N,1); % mean = x
        v_tau = 100*tau*ones(size(y));
        v_tau(k) = tau;
        [~,prio_xY,post_xY] = forward_recursion(P, p_x(:,1),v_tau,y,N);
        total_post_xY = backward_recursion(P,prio_xY,post_xY,N);
        value_dont_clean = -5000*sum(total_post_xY(2,:));
        mean_PoV = mean_PoV + max(value_clean_all,value_dont_clean)/mc_samples;
    end
    PoV(k) = mean_PoV;
end
VoI_mc = PoV - PV;
maximum = max(VoI_mc)
max_location = find(VoI_mc==maximum)
%%
figure(3); hold off;
plot(1:N,VoI_mc)
xlabel('railroad section $k$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('VOI$(k)$', 'interpreter', 'latex', 'FontSize', 16);
%legend('$p(x_i = 1|\mathbf{y}_D)$','$\mathbf{y}_D$', 'Location', 'best', 'interpreter', 'latex', 'FontSize', 14);
grid

%% Numerical integration
% Much faster and more accurate than MC, but currently 10 to low values

%p(y_k | x_k = x) ~ N(x,tau); N(mean, standard_deviation) p(y_k) = sum_{x
%in {0,1}} ( p(y_k | x_k = x)*p(x_k = x) ) p(y_k) ~ N(0, tau)*p(x_k=0) +
%N(1, tau)*p(x_k=1)   where p(x_k) >= 0 p(y_k) ~ N(0, tau*p(x_k=0)) +
%N(p(x_k=1), tau*p(x_k=1)) p(y_k) ~ N(p(x_k=1),sqrt((tau*p(x_k=0))^2 +
%(tau*p(x_k=1))^2)) p(y_k) ~ N(p(x_k=1),tau*sqrt(p(x_k=0)^2 + p(x_k=1)^2))

Nsamples = 101;
spread = 4; %standard devitaions 
mu = p_x(2,:)';
std = tau*sqrt(sum(p_x.^2,1))';
step_size = 2*spread*std/(Nsamples-1);
Ysamples = zeros(N,Nsamples);
for k = 1:N
    Ysamples(k,:) = (mu(k)-2*spread*std(k):step_size(k):mu(k)) + spread*std(k);
end
% discretized pdf: dpdf_y(k,i) = p(y_k = Ysamples(k,i))
dpdf_y = normpdf(Ysamples, mu , std);

figure
hold on
plot(Ysamples(1,:),dpdf_y(1,:),'*')
plot(Ysamples(N,:),dpdf_y(N,:),'*')
grid
legend("p(y_1)","p(y_N)")

PoV = -ones(N,1);
for k = 1:k
    PoV_ni = 0;
    for i = 1:Nsamples
        y = normrnd(x,tau,N,1);
        v_tau = 100*tau*ones(size(y));
        v_tau(k) = tau;
        y(k) = Ysamples(k,i);
        [~,prio_xY,post_xY] = forward_recursion(P, p_x(:,1),v_tau,y,N);
        total_post_xY = backward_recursion(P,prio_xY,post_xY,N);
        value_dont_clean = -5000*sum(total_post_xY(2,:));
        integrand_PoV = max(value_clean_all,value_dont_clean)*dpdf_y(k,i);
        PoV_ni = PoV_ni + integrand_PoV * step_size(k);
    end
    PoV(k) = PoV_ni;
end
VoI_ni = PoV - PV;    
figure
plot(VoI_ni)
xlabel('locations')
ylabel('VoI')
legend("Numerical integration")
fig4 = figure;
hold on
plot(1:N,VoI_mc)
plot(VoI_ni)
legend("Monte carlo","Numerical integration")
grid
xlabel('locations')
ylabel('VoI')
title_pic = 'VoI';
% %save_plots(title_pic,directory_d,fig4);

%% task e
big_pov = zeros(N,N);

% diagonal
for k=1:N
    mean_PoV = 0;
    for j=1:mc_samples
        x_new = draw_from_HMM(P,p_x(:,1),N);
        y = normrnd(x_new,tau/sqrt(2),N,1); % mean = x
        v_tau = 100*tau*ones(size(y));
        v_tau(k) = tau;
        [~,prio_xY,post_xY] = forward_recursion(P, p_x(:,1),v_tau,y,N);
        total_post_xY = backward_recursion(P,prio_xY,post_xY,N);
        value_dont_clean = -5000*sum(total_post_xY(2,:));
        mean_PoV = mean_PoV + max(value_clean_all,value_dont_clean)/mc_samples;
    end
    big_pov(k,k) = mean_PoV;
end

% off-diagonal
for k=1:N
    for j=k:N
        mean_PoV = 0;
        if(j==k)
            continue;
        end
        for i=1:mc_samples
            x_new = draw_from_HMM(P,p_x(:,1),N);
            y = normrnd(x_new,tau,N,1); % mean = x
            v_tau = 100*tau*ones(size(y));
            v_tau(k) = tau;
            v_tau(j) = tau;
            [~,prio_xY,post_xY] = forward_recursion(P, p_x(:,1),v_tau,y,N);
            total_post_xY = backward_recursion(P,prio_xY,post_xY,N);
            value_dont_clean = -5000*sum(total_post_xY(2,:));
            mean_PoV = mean_PoV + max(value_clean_all,value_dont_clean)/mc_samples;
        end
        big_pov(k,j) = mean_PoV;
        big_pov(j,k) = big_pov(k,j);
    end
end

big_voi = big_pov - PV;
maximum = max(max(big_voi));
[max_x,max_y]=find(big_voi==maximum)
%%
fig5 = figure;
imagesc(big_voi)
set(gca,'YDir','normal')
xlabel('Sensor Placement 1')
ylabel('Sensor Placement 2')
colorbar
title_pic = 'VoI_cross';
% save_plots(title_pic,directory_e,fig5);

%% task f
mc_samples = 10000;
PoV_f = zeros(mc_samples,1);
tau_new = 1;

for k=1:mc_samples
    x_new_f = draw_from_HMM(P,p_x(:,1),N);
    y = normrnd(x_new_f,tau_new,N,1);
    v_tau = ones(size(y));
    [~,prio_xY,post_xY] = forward_recursion(P, p_x(:,1),v_tau,y,N);
    total_post_xY = backward_recursion(P,prio_xY,post_xY,N);
    value_dont_clean = -5000*sum(total_post_xY(2,:));
    %mean_PoV = mean_PoV + max(value_clean_all,value_dont_clean)/mc_samples;
    PoV_f(k) = max(value_clean_all,value_dont_clean);
end

VoI_f = PoV_f - PV;
%%
fig6 = figure;
%plot(1:mc_samples,VoI_f,'*')
histogram(VoI_f)
set(gca,'YScale','log')
xlabel('VOI')
ylabel('# MC samples')
grid
title_pic = 'VoI_f';
% save_plots(title_pic,directory_f,fig6);





