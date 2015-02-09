figure(1);hold all;
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
figure(2); hold all;