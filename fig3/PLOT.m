figure
p1=plot(((A_max-A_min)/(N_A-1)*(0:N_A-1) + A_min)/lambda,SE_average,'s-m','linewidth',2);
hold on
p2=plot(((A_max-A_min)/(N_A-1)*(0:N_A-1) + A_min)/lambda,SE_base,'o-b','linewidth',2);
hold on 
p3=plot(((A_max-A_min)/(N_A-1)*(0:N_A-1) + A_min)/lambda,SE_base,'d-r','linewidth',2);
hold on
legend([p1 p2 p3],'ROMA-MIMOS ','MA-MIMOS','FAS','Interpreter','latex' )

xlabel('Normalized region size ($\lambda$)','Interpreter','latex')
ylabel('Spectral Efficiency (bits/s/Hz)','Interpreter','latex')
grid on