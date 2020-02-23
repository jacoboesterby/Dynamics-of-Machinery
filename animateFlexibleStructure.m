function animateFlexibleStructure(x,phi,t_int)
close all
m1 = 2.294;
m2 = 1.941;
m3 = 1.943;
m4 = 2.732;
m5 = 0.774;
m6 = 0.774;
mbeam = 0.3;
rho1 = 0.278;
mbeam2 = 0.421;
rho2 = 0.401;
macc = 0.073;

E = 2*10^11;
b1 = 0.0291;
h1 = 0.0011;
I1 = b1*h1^3/12;
L1 = 0.193;
L2 = 0.292;
L3 = 0.220;
L4 = 0.271;
L5 = 0.428;
L6 = 0.428;

b2 = 0.0350;
h2 = 0.0015;
I2 = b2*h2^3/12;


Sigma = 0.7341;
Phi = @(xi,L) cosh(1.87510407/L*(L-xi))+cos(1.87510407/L*(L-xi))-Sigma*(sinh(1.87510407/L*(L-xi))+sin(1.87510407/L*(L-xi)));
%Define system of equations 
As = diag(ones(1,6).*2);
g = 9.81;
%%Plotting stuff 
fps = 60;
time = linspace(0,t_int(end),t_int(end)*fps);
%Create .avi file
f = figure('visible','off');
v = VideoWriter('Test5.avi');
open(v);
hold on
Beam = @(L) linspace(0,L,20);

sol =@(s) interp1(t_int.',x.',s);
ang =@(s) interp1(t_int.',phi.',s);
w = 0.2;
h = 0.04;
rhub = 0.05;


for i=time
    q = linsolve(As,sol(i).');
    %Beam and Mass 1
    x1 = q(1)*Phi(Beam(L1),L1);
    y1 = Beam(L1);
    plot(x1+w,y1,'linewidth',3,'color','k')
    plot(x1-w,y1,'linewidth',3,'color','k')
    
    rectangle('Position',[x1(end)-w,y1(end),w*2,h],'FaceColor','k',...
    'LineWidth',3)

    %Beam and Mass 2
    x2 = q(2)*Phi(Beam(L2),L2)+x1(end);
    y2 = Beam(L2)+y1(end)+h;
    plot(x2+w,y2,'linewidth',3,'color','k')
    plot(x2-w,y2,'linewidth',3,'color','k')
    
    rectangle('Position',[x2(end)-w,y2(end),w*2,h],'FaceColor','k',...
    'LineWidth',3)
    
    %Beam and Mass 3
    x3 = q(3)*Phi(Beam(L3),L3)+x2(end);
    y3 = Beam(L3)+y2(end)+h;
    plot(x3+w,y3,'linewidth',3,'color','k')
    plot(x3-w,y3,'linewidth',3,'color','k')
    
    rectangle('Position',[x3(end)-w,y3(end),w*2,h],'FaceColor','k',...
    'LineWidth',3)

    %Beam and Mass 4
    x4 = q(4)*Phi(Beam(L4),L4)+x3(end);
    y4 = Beam(L4)+y3(end)+h;
    plot(x4+w,y4,'linewidth',3,'color','k')
    plot(x4-w,y4,'linewidth',3,'color','k')
    
    rectangle('Position',[x4(end)-w,y4(end),w*2,h],'FaceColor','k',...
    'LineWidth',3)
    
    %Blades
    % 1
    x5 = (q(5)*Phi(Beam(L5),L5))*cos(ang(i)) +  (Beam(L5))*(-sin(ang(i))) + x4(end) - sin(ang(i))*rhub;
    y5 = (q(5)*Phi(Beam(L5),L5))*sin(ang(i)) +  (Beam(L5))*cos(ang(i)) + y4(end) + h + cos(ang(i))*rhub;
    plot(x5,y5,'linewidth',3,'color','m')
    plot(x5,y5,'linewidth',3,'color','m')
    
    % 2
    x6 = (-q(6)*Phi(Beam(L6),L6))*cos(ang(i)-pi) +  (Beam(L5))*(-sin(ang(i)-pi)) + x4(end) + sin(ang(i))*rhub;
    y6 = (-q(6)*Phi(Beam(L6),L6))*sin(ang(i)-pi) +  (Beam(L5))*cos(ang(i)-pi) + y4(end) + h - cos(ang(i))*rhub;
    plot(x6,y6,'linewidth',3,'color','m')
    plot(x6,y6,'linewidth',3,'color','m')
    
    filledCircle([x4(end),y4(end)+h],rhub,100,'g');
    
    xlim([-w*3,w*3]);
    ylim([0,L1+L2+L3+L4+L6].*1.1)
    frame = getframe(gcf);
    writeVideo(v,frame)
    fprintf('Progress: %.2f\n',i/time(end)*100);
    cla;
end
close(v);