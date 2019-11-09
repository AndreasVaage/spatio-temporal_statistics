close all
clear

%% task a
data = 50;
marg_p = zeros(data,1);
marg_p(1) = 0.01;

for i=2:data
    marg_p(i)= marg_p(i-1) + 0.05*(1-marg_p(i-1));
end

fig1 = figure;
plot(1:data,marg_p)
grid
xlabel('locations')
ylabel('marg prob')

%% task c
tau = 0.3;
x = zeros(data,1);
%x(1) = binornd(1,marg_p(1));
if(rand>=0.99)
    x(1) = 1;
else
    x(1) = 0;
end

for i=2:data
   if(x(i-1)==1)
       x(i) = 1;
   else
       if(rand<0.05)
           x(i) = 1;
       else
           x(i) = 0;
       end
   end
end

y = normrnd(x,tau,data,1); % same as mu + sigma*randn(data,1)
y(20) = 0.2;
y(30) = 0.7;

% fig2 = figure;
% plot(1:data,y);
% grid

[~,prio_xY,post_xY] = forward_recursion(0.95,tau,y,data);
margp_xY = backward_recursion(0.95,prio_xY,post_xY,data);

fig3 = figure;
plot(1:data,margp_xY(2,:));
grid
xlabel('locations')
ylabel('posterior prob for x_i = 1') 
% Why this is that shitty??

%% task d
PoV = [];
mc_samples = 10000;
clean_all = -100000;
PV = -100000;

for i=1:data
    pov = 0; 
    for j=1:mc_samples
        y = normrnd(x,100*tau,data,1); % what mean?
        y(i) = normrnd(binornd(1,marg_p(i)),tau);
        [~,prio_xY,post_xY] = forward_recursion_d(0.95,tau,y,data,i);
        prob = backward_recursion(0.95,prio_xY,post_xY,data);
        dont_clean = -5000*sum(prob(2,:));
        pov = pov + max(clean_all,dont_clean)/j;
    end
    PoV(i) = pov;
end
%PoV = PoV(:,any(PoV));
VoI = PoV - PV

% fig4 = figure;
% plot(1:data,VoI)
% grid
% xlabel('locations')
% ylabel('VoI') 


%% task e
big_pov = zeros(data,data);

% diagonal
for i=1:data
    for j=1:mc_samples
        y = normrnd(0.5,100*tau,data,1);
        x_new = zeros(data,1);
        x_new(i) = binornd(1,marg_p(i));
        y(i) = normrnd(x_new(i),tau/sqrt(2));
        [~,prio_xY,post_xY] = forward_recursion_d(0.95,tau,y,data,i);
        prob = backward_recursion(0.95,prio_xY,post_xY,data);
        dont_clean = -5000*sum(prob(2,:));
        pov = pov + max(clean_all,dont_clean)/j;
    end
    big_pov(i,i) = pov;
end

% off-diagonal
for i=10:data
    for j=i:data
        if(j==i)
            break;
        end
        for k=1:mc_samples
            y = normrnd(0.5,100*tau,data,1);
            x_new = zeros(data,1);
            x_new(i) = binornd(1,marg_p(i));
            for l=i+1:j
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
            y(i) = normrnd(x_new(i),tau);
            y(j) = normrnd(x_new(i),tau);
            [~,prio_xY,post_xY] = forward_recursion_e(0.95,tau,y,data,i,j);
            prob = backward_recursion(0.95,prio_xY,post_xY,data);
            dont_clean = -5000*sum(prob(2,:));
            pov = pov + max(clean_all,dont_clean)/k;
        end
    end
    big_pov(i,j) = pov;
    big_pov(j,i) = big_pov(i,j);
end

big_voi = big_pov - PV;
maximum = max(max(big_voi));
[max_x,max_y]=find(big_voi==maximum)

% plot matrix.
% plot diagonal: 2 sensors in same loc.
% plot 20,30

%% task f
PoV_f = zeros(mc_samples,1);
tau_new = 1;

for i=1:mc_samples
    x_new_f = zeros(data,1);
    if(rand>=0.99)
        x(1) = 1;
    else
        x(1) = 0;
    end

    for i=2:data
       if(x(i-1)==1)
           x(i) = 1;
       else
           if(rand<0.05)
               x(i) = 1;
           else
               x(i) = 0;
           end
       end
    end
    y = normrnd(x_new_f,tau_new,data,1);
    [~,prio_xY,post_xY] = forward_recursion_f(0.95,tau_new,y,data);
    prob = backward_recursion(0.95,prio_xY,post_xY,data);
    dont_clean = -5000*sum(prob(2,:));
    PoV_f(i) = max(clean_all,dont_clean);
end

VoI_f = PoV_f - PV;

% plot some of the new VoI: should we clean the road upfront or not?






