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
        y(i) = normrnd(x(i),tau);
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
% plot(1:36,VoI)
% grid
% xlabel('locations')
% ylabel('VoI') 


%% task e












