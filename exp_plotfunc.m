

function exp_plotfunc(x,y,t,F,H1y,H2y,fase1y,fase2y,Cy,idx)


%_____________________________________________________
% Graphics 1 - Analysis in the time domain
figure(1)

subplot(2,1,1); plot(t,y,'k-'); 
grid
title('(a) Acceleration in Time Domain')
ylabel('(a) acc [m/s^2]')
xlabel('time [s]')
axis([0 max(t) -5 5])
subplot(2,1,2); plot(t,x,'k-');
grid
title('(b) Force Signal in Time Domain')
ylabel('(b) force [N]')
xlabel('time [s]')
axis([0 max(t) -5 5])


%_____________________________________________________
% Graphics 2 - Analysis in the frequency domain
figure(2)

subplot(3,1,1), plot(F(idx),Cy(idx),'k-','LineWidth',1.5); 
axis([0 max(F(idx)) 0 1.1])
set(gca,'FontAngle','oblique','FontSize',14)
ylabel('Coherence ')
grid

subplot(3,1,2), plot(F(idx),abs(H1y(idx)),'b-',F(idx),abs(H2y(idx)),'r-','LineWidth',1.5); 
axis([0 max(F(idx)) 0  1.1*max(abs(H2y(idx)))])
grid
set(gca,'FontAngle','oblique','FontSize',14)
ylabel('FRF [(m/s^2)/N] ')
legend('H1','H2')
subplot(3,1,3), plot(F(idx),fase1y(idx)*180/pi,'b-',F(idx),fase2y(idx)*180/pi,'r-','LineWidth',1.5); 
axis([0 max(F(idx)) 0 200])
grid
set(gca,'FontAngle','oblique','FontSize',14,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
set(gca,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
axis([0 max(F(idx)) 0 200])
xlabel('Frequency [Hz]')
ylabel('Phase [^o]')
legend('H1','H2')

%_____________________________________________________
% Graphics 3 - Analysis in the frequency domain - H2
figure(3)

subplot(2,1,1), semilogy(F(idx),abs(H2y(idx)),'r-','LineWidth',1.5); 
axis([0.7 max(F(idx)) 0  1.1*max(abs(H2y(idx)))])
grid
set(gca,'FontAngle','oblique','FontSize',14)
ylabel('FRF [(m/s^2)/N] ')
subplot(2,1,2), plot(F(idx),fase2y(idx)*180/pi,'r-','LineWidth',1.5); 
grid
set(gca,'FontAngle','oblique','FontSize',14,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
set(gca,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
axis([0.7 max(F(idx)) 0 200])
xlabel('Frequency [Hz]')
ylabel('Phase [^o]')


