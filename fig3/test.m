y = [5.20941, 4.83458, 4.76554, 4.79881, 4.76085];
figure(1);
h3=plot(((A_max-A_min)/(N_A-1)*(0:N_A-1) + A_min)/lambda, y,'-','linewidth',1.5);
hold on;
