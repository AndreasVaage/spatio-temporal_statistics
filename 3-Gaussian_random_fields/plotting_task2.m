function plotting_task2(theta, alpha, LogLikelihood, alpha_true, theta_true)

figure(3)
hold off;
subplot(5,1,1); hold on;
title('\textbf{10 iterations}', 'Interpreter', 'latex', 'FontSize', 15);
yline(alpha_true,'-.');
plot(alpha)
legend("$\hat{\alpha}$","$\alpha$",'interpreter', 'latex', 'FontSize', 15)
hold off;
subplot(5,1,2); hold on;
yline(theta_true(1),'-.');
plot(theta(1,:))
legend("$\hat{\sigma}^2$","$\sigma^2$",'interpreter', 'latex', 'FontSize', 15)
hold off;
subplot(5,1,3); hold on;
yline(theta_true(2),'-.');
plot(theta(2,:))
legend("$\hat{\eta}$","$\eta$",'interpreter', 'latex', 'FontSize', 15)
subplot(5,1,4); hold on;
yline(theta_true(3),'-.');
plot(theta(3,:))
legend("$\hat{\tau}^2$","$\tau^2$",'interpreter', 'latex', 'FontSize', 15)
subplot(5,1,5); hold on;
plot((LogLikelihood))
legend("$l(\mathbf{Y} | \hat\theta,\hat\alpha)$",'interpreter', 'latex', 'FontSize', 15)
hold off;
end