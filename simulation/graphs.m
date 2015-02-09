figure();hold all;
plot([0 .47],[0 .62],'kx:','markersize',10);
plot([.47 1],[.62 0],'k:');
% rint
plot([0 1],[0 0],'kx-','markersize',10);
% rtil
plot([0 .3],[0 0],'r--','linewidth',2);
text(.15,-.1,'$\tilde{r}$','interpreter','latex');
% Rtil
plot([.3 .47],[0 .62],'b-','linewidth',2);
text(.42,.34,'$\tilde{R}$','interpreter','latex');
% theta
ang([.3 0],.08,[0 1.2],'k');
text(.37,.08,'$\theta$','interpreter','latex');

xlim([-1 2]); ylim([-1 1]);
set(gca,'visible','off');
saveas(gcf,'LocationDensityExplain.png')
%%
data = dir('*CrashRadius.mat');
data = {data.name};
figure(); hold all; grid on;
load(data{1})
plot(R/1e3,PR);
load(data{2})
plot(R/1e3,PR);
load(data{3})
plot(R/1e3,PR);
legend('Airbus 380','B737-900ER','G280','location','best')
xlim([0 300]); xlabel('Crash Radius [km]'); ylabel('Probability')
title('Estimated Crash Radius Distribution')
saveas(gcf,'CrashRadiusExplain.png');
%%
data = dir('nosearchEsc1*.mat');
data = {data.name};
figure(); hold all; grid on;
load(data{1})
plot(tVec,qt);
load(data{2})
plot(tVec,qt);
legend('5x5[km^2] Grid','20x15[km^2] Grid','location','northwest')
ylim([0 1]); xlabel('Tiem [hr]'); ylabel('Probability')
title('Probability of Aircraft Debris Escaping the Search Domain')
saveas(gcf,'NoSearchEscapeExplain.png');