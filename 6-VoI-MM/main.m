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
figure
hold on
plot(p_x(1,:))
plot(p_x(2,:))
fig1 = figure;
plot(1:N,marg_p)
grid
xlabel('locations')
ylabel('marg prob')
title_pic = 'marg_prob';
%save_plots(title_pic,directory_a,fig1);
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
% fig2 = figure;
% plot(1:data,y);
% grid

[~,prio_xY,post_xY] = forward_recursion(P,p_x(:,1),v_tau,y,N);
margp_xY = backward_recursion(P,prio_xY,post_xY,N);

fig3 = figure;
plot(1:N,margp_xY(2,:));
grid
xlabel('locations')
ylabel('posterior prob for x_i = 1')
title_pic = 'post_prob';
%save_plots(title_pic,directory_c,fig3);

%% task d
value_clean_all = -100000;
PV = -100000;

%% Monte carlo
PoV = -ones(N,1);
mc_samples = 1000;
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
    
fig4 = figure;
hold on
plot(1:N,VoI_mc)
grid
xlabel('locations')
ylabel('VoI')
title_pic = 'VoI';
%save_plots(title_pic,directory_d,fig4);


%% task e
big_pov = zeros(N,N);

% diagonal
for k=1:N
    for j=1:mc_samples
        y = normrnd(0.5,100*tau,N,1);
        x_new = zeros(N,1);
        x_new(k) = binornd(1,marg_p(k));
        y(k) = normrnd(x_new(k),tau/sqrt(2));
        [~,prio_xY,post_xY] = forward_recursion_d(0.95,tau,y,N,k);
        prob = backward_recursion(0.95,prio_xY,post_xY,N);
        dont_clean = -5000*sum(prob(2,:));
        pov = pov + max(value_clean_all,dont_clean)/j;
    end
    big_pov(k,k) = pov;
end

% off-diagonal
for k=10:N
    for j=k:N
        if(j==k)
            break;
        end
        for k=1:mc_samples
            y = normrnd(0.5,100*tau,N,1);
            x_new = zeros(N,1);
            x_new(k) = binornd(1,marg_p(k));
            for l=k+1:j
                if(x_new(l-1) == 1)
                    x_new(l) = 1;
                else
                    if(rand(1) < 0.05)
                        x(l) = 1;
                    else
                        x(l) = 0;
                    end
                end
            end
            y(k) = normrnd(x_new(k),tau);
            y(j) = normrnd(x_new(k),tau);
            [~,prio_xY,post_xY] = forward_recursion_e(0.95,tau,y,N,k,j);
            prob = backward_recursion(0.95,prio_xY,post_xY,N);
            dont_clean = -5000*sum(prob(2,:));
            pov = pov + max(value_clean_all,dont_clean)/k;
        end
    end
    big_pov(k,j) = pov;
    big_pov(j,k) = big_pov(k,j);
end

big_voi = big_pov - PV;
maximum = max(max(big_voi));
[max_x,max_y]=find(big_voi==maximum)

% plot matrix.
% plot diagonal: 2 sensors in same loc.
% plot 20,30
% save plots in directory_e.

%% task f
PoV_f = zeros(mc_samples,1);
tau_new = 1;

for k=1:mc_samples
    x_new_f = zeros(N,1);
    if(rand>=0.99)
        x(1) = 1;
    else
        x(1) = 0;
    end

    for k=2:N
       if(x(k-1)==1)
           x(k) = 1;
       else
           if(rand<0.05)
               x(k) = 1;
           else
               x(k) = 0;
           end
       end
    end
    y = normrnd(x_new_f,tau_new,N,1);
    [~,prio_xY,post_xY] = forward_recursion_f(0.95,tau_new,y,N);
    prob = backward_recursion(0.95,prio_xY,post_xY,N);
    dont_clean = -5000*sum(prob(2,:));
    PoV_f(k) = max(value_clean_all,dont_clean);
end

VoI_f = PoV_f - PV;

% plot some of the new VoI: should we clean the road upfront or not?
% save plots in directory_f.





