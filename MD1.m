clear all
syms E I1 I2 L L1 L2 L3 L4 L5 L6 m1 m2 m3 m4 m5 m6 x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) t g theta(t) R1(t) R2(t) R3(t) R4(t) R5(t) R6(t) u1 u2 u3 u4 u5 u6 phi

eq1 = sym(zeros(3,1));
eq2 = sym(zeros(3,1));
eq3 = sym(zeros(3,1));
eq4 = sym(zeros(3,1));
eq5 = sym(zeros(3,1));
eq6 = sym(zeros(3,1));

%Position vectors
iroa = [x1(t),L1,0].';
b1rab = [x2(t),L2,0].';
b2rbc = [x3(t),L3,0].';
b3rcd = [x4(t),L4,0].';
b5rde = [x5(t),L5,0].';
b5rdf = [x6(t),-L6,0].';

%Transformation matrices:
T_theta = [cos(theta), sin(theta), 0;
    -sin(theta), cos(theta),0;
    0,           0,         1];
T = [1,0,0;
    0,1,0;
    0,0,1];

omega = [0,0,diff(theta,t)].';

%Velocity vectors:
b1va = diff(iroa,t);
b2vb = b1va + diff(b1rab,t);
b3vc = b2vb + diff(b2rbc,t);
b4vd = b3vc + diff(b3rcd,t);

%Velocity of blade tips
b5ve = T_theta*b4vd + cross(omega,b5rde) + diff(b5rde,t);
b5vf = T_theta*b4vd + cross(omega,b5rdf) + diff(b5rdf,t);



%Acceleration vectors:
%Note that omega = [0,0,0] for all masses m1 m2 m3 m4, such that all cross
%products including omega can be neglected
b1Aa = diff(iroa,t,t);
b2Ab = b1Aa + diff(b1rab,t,t);
b3Ac = b2Ab + diff(b2rbc,t,t);
b4Ad = b3Ac + diff(b3rcd,t,t);

%Accelerations for blade tips:
b5Ae = T_theta*b4Ad + cross(diff(omega,t),b5rde) + cross(omega,cross(omega,b5rde)) + 2 * cross(omega,diff(b5rde,t)) + diff(b5rde,t,t);
b5Af = T_theta*b4Ad + cross(diff(omega,t),b5rdf) + cross(omega,cross(omega,b5rdf)) + 2 * cross(omega,diff(b5rdf,t)) + diff(b5rdf,t,t);



BeamStiff_blade = @(L,I) 3*E*I/(L^3);
BeamStiff_beam = @(L,I) 12*E*I/(L^3);

IPA = [0, -m1*g, 0].';
ITA1 = [-2*BeamStiff_beam(L1,I1)*x1(t), 0, 0].';
IRA1 = [0,R1,0].';

ITA2 = [2*BeamStiff_beam(L2,I1)*(x2(t)-x1(t)),0,0].';
IRA2 = [0,-R2,0].';

IPB = [0, -m2*g, 0].';
ITB2 = -ITA2;
IRB2 = -IRA2;

ITB3 = [(x3(t)-x2(t))*BeamStiff_beam(L3,I1)*2,0,0].';
IRB3 = [0,-R3,0].';

IPC = [0,-m3*g,0].';
ITC3 = -ITB3;
IRC3 = -IRB3;
ITC4 = [(x4(t)-x3(t))*BeamStiff_beam(L4,I1)*2,0,0].';
IRC4 = [0,-R4,0].';

IPD = [0,-m4*g,0].';
ITD4 = -ITC4;
IRD4 = -IRC4;
%


IPE = T_theta*[0,-m5*g,0].';
ITE = [-BeamStiff_blade(L5,I2)*x5(t),0,0].';
IRE = [0,R5,0].';



IPF = T_theta*[0,-m6*g,0].';
ITF = [-BeamStiff_blade(L6,I2)*x6(t),0,0].';
IRF = [0,R6,0].';


IRD5 = T_theta.'*(-(IRE+IRF));
ITD5 = T_theta.'*(-(ITE+ITF));

ITEF5 = ITE+ITF;

%ITEF5 =

eq1(:) = m1.*(b1Aa) - ( IPA + IRA1 + IRA2 + ITA1 + ITA2);
eq2(:) = m2.*(b2Ab) - ( IPB + IRB2 + IRB3 + ITB2 + ITB3);
eq3(:) = m3.*(b3Ac) - ( IPC + IRC3 + IRC4 + ITC3 + ITC4);
eq4(:) = m4.*(b4Ad) - ( IPD + IRD4 + IRD5 + ITD4 + ITD5);
eq5(:) = m5.*(b5Ae) - ( IPE + IRE + ITE);
eq6(:) = m6.*(b5Af) - ( IPF + IRF + ITF);

% eq1(:) = m1.*(b1Aa) - ( IPA+ITA1+ITA2);
% eq2(:) = m2.*(b2Ab) - ( IPB + ITB2 + ITB3);
% eq3(:) = m3.*(b3Ac) - ( IPC + ITC3 + ITC4);
% eq4(:) = m4.*(b4Ad) - ( IPD + ITD4 + ITD5);
% eq5(:) = m5.*(b5Ae) - ( IPE + ITE - ITD5);
% eq6(:) = m6.*(b5Af) - ( IPF + ITF - ITD5);


%Numerical values:
m1 = 2.294;
m2 = 1.941;
m3 = 1.943;
m4 = 2.732;
m5 = 0.774;
m6 = 0.774+0.02;
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

BeamStiff_blade = @(L,I) 3*E*I/(L^3);
BeamStiff_beam = @(L,I) 12*E*I/(L^3);

b2 = 0.0350;
h2 = 0.0015;
I2 = b2*h2^3/12;

%theta = pi/4;

g = 9.81;


%Analytical eigen value analysis
%Mass matrix
M = [m1,0,0,0,0,0;
     0,m2,0,0,0,0;
     0,0,m3,0,0,0;
     0,0,0,m4,0,0;
     0,0,0,0,m5,0;
     0,0,0,0,0,m6];
 
k1 = BeamStiff_beam(L1,I1);
k2 = BeamStiff_beam(L2,I1);
k3 = BeamStiff_beam(L3,I1);
k4 = BeamStiff_beam(L4,I1);
k5 = BeamStiff_blade(L5,I2);
k6 = BeamStiff_blade(L6,I2);


K = [2*k1+2*k2  -2*k2      0         0         0   0;
     -2*k2 2*k2+2*k3  -2*k3     0         0   0;
     0     -2*k3      2*k3+2*k4 -2*k4     0   0;
     0     0          -2*k4     2*k4+k5+k6 -k5 -k6;
     0     0          0         -k5       k5  0;
     0     0          0         -k6       0   k6];
nulmat = zeros(size(M));


A = [M,nulmat;
     nulmat,M];
B = [nulmat,K;
     -M,nulmat]; 
 
 
 %%% Bode plot
% 
% Sys = linsolve(-B,A);
% U = [ones(6,1);zeros(6,1)]
% C  = diag(ones(1,12));
% C(7:end,:) = 0;
% D  = zeros(12,1);
% Sys = ss(Sys,U,C,D)
% 
% bode(Sys)

 
 %%%
 
 
[u,lambda] = eig(-B,A);

lambda = diag(lambda);

[lambda,ind] = sort(lambda);
omega0 = lambda
f0 = lambda./(2*pi);
u = u(:,ind);
Phi = u(1:6,1:2:end);
%PlotModeShapes(Phi,f0(2:2:end));

%% Damping
xi =[0.0009,0.0007,0.0009,0.0013,0.001,0.0013].';
xi = [0.0024,0.0018,0.0020,0.0031,0.0028,0.0040].';




%Solve for all damping ratios: 

for i=1:6
    temp1(i) = 1/(2*omega0(i*2-1));
    temp2(i) = omega0(i*2-1)/2;
end
C = [temp1.',temp2.'];
res = linsolve(C,xi);
res = (C.'*C)^(-1)*C.'*xi;
alpha = abs(res(1));
beta = abs(res(2));


syms ddx ddx1 ddx2 ddx3 ddx4 ddx5 ddx6 dx x dx1 dx2 dx3 dx4 dx5 dx6 xx1 xx2 xx3 xx4 xx5 xx6
x_ddot = [ddx1,ddx2,ddx3,ddx4,ddx5,ddx6].';
x_dot = [dx1,dx2,dx3,dx4,dx5,dx6].';
x = [xx1,xx2,xx3,xx4,xx5,xx6].';

D1 = alpha.*M;
D2 = beta.'*K;
D3 = alpha.*M+beta.*K;

A = [M,D3;
     nulmat,M];
B = [nulmat,K;
     -M,nulmat]; 
 
z = sym(zeros(12,1));
z_dot = sym(zeros(12,1));
z(:) = [x_dot;x];
z_dot(:) = [x_ddot;x_dot];

[u,lambda] = eig(-B,A);

lambda = diag(lambda);

[lambda,ind] = sort(lambda);
omegad = lambda;
fd = omegad./(2*pi);
u = u(:,ind);
Phi = u(1:6,1:2:end);

for i =1:12
    xi(i) = -real(omegad(i))/sqrt((real(omegad(i)))^2+(imag(omegad(i)))^2);
end

%PlotModeShapes(Phi,fd(2:2:end));

Sys = inv(A)*(-B);

U = [zeros(6,1);ones(1,1);zeros(5,1)];

C =  zeros(12,12);
for i=1:12
    C(i,i) = 1;
end
D = zeros(12,1);

sys = ss(Sys,U,C,D);
w = linspace(0,8*2*pi,1000);
[mag,phase,fr] = bode(sys,w);
amp = squeeze(mag(7,1,:));
ph = squeeze(phase(7,1,:));
figure
hold on
subplot(2,1,1)
grid on
semilogy(fr/(2*pi),amp)
%subplot(2,1,2)
%plot(fr/(2*pi),ph+90);

srate = 50;
sf = 1/srate;

f1 = struct;
load('force1.txt');
f1.signals.values = force1;
f1.time = linspace(sf,length(force1)*sf,length(force1));

f2 = struct;
load('force2.txt');
f2.signals.values = force2;
f2.time = linspace(sf,length(force2)*sf,length(force2));

f3 = struct;
load('force3.txt');
f3.signals.values = force3;
f3.time = linspace(sf,length(force3)*sf,length(force3));

f4 = struct;
load('force4.txt');
f4.signals.values = force4;
f4.time = linspace(sf,length(force4)*sf,length(force4));

f5 = struct;
load('force5.txt');
f5.signals.values = force5;
f5.time = linspace(sf,length(force5)*sf,length(force5));

f6 = struct;
load('force6.txt');
f6.signals.values = force6;
f6.time = linspace(sf,length(force6)*sf,length(force6));


sol = sim('SimulinkModelDamped',length(force1)*sf);
plot(sol.tout,sol.x2sim,'linewidth',3)
[sol1,time] = step(sys,160);
hold on
plot(time,sol1(:,2),'linestyle','--','linewidth',3)

% [Pxx,F,Pxxc] = pwelch(sol1(:,2),[],[],[],8*2*pi,'twosided')
% figure
% plot(F,Pxx)
% 
% pwelch(sol1(:,2))

load('acceleration1.txt');
load('force1.txt')
GenerateFRF(acceleration1,force1);
GenerateFRF(sol.xdd1sim(1:length(force1)),force1);
%pwelch, cpsd

%%%%%%%%%%%




%%

eq1 = subs(eq1);
eq2 = subs(eq2);
eq3 = subs(eq3);
eq4 = subs(eq4);
eq5 = subs(eq5);
eq6 = subs(eq6);

eqsys = [eq1;eq2;eq3;eq4;eq5;eq6];
syms ddx1 ddx2 ddx3 ddx4 ddx5 ddx6 RA RB RC RD RE RF
dderivatives = [ddx1 ddx2 ddx3 ddx4 ddx5 ddx6];
forces = [RA RB RC RD RE RF];
eqsys = subs(eqsys,[diff(x1,t,t),diff(x2,t,t),diff(x3,t,t),diff(x4,t,t),diff(x5,t,t),diff(x6,t,t),R1,R2,R3,R4,R5,R6],[dderivatives,forces]);
dx = [ddx1,ddx2,ddx3,ddx4,ddx5,ddx6,RA,RB,RC,RD,RE,RF];



index = [];

for i=1:size(eqsys,1)
    if mod(i,3) ==0
        
    else
        index(end+1) = i;
    end
end

keyboard



[A,b] = equationsToMatrix(eqsys(index),dx.');
syms ddphi

%b = subs(b,{x1,x2,x3,x4,x5,x6},{xx1,xx2,xx3,xx4,xx5,xx6});
%[C,d] = equationsToMatrix(b,[xx1,xx2,xx3,xx4,xx5,xx6]);


% numerical solution
acc = linsolve(A,b);
acc = acc(1:6);


u = sym(zeros(6,1));
u = [u1 u2 u3 u4 u5 u6];
ddphi = 0;
acc = subs(acc,[x1,x2,x3,x4,x5,x6,theta,diff(theta,t,t)],[u(1),u(2),u(3),u(4),u(5),u(6),phi,ddphi]);


acc = matlabFunction(acc);
%%
clear t_int theta
% time step
deltaT = 0.00001;

% number of integration points
n_int = 500000;

%Initial conditions on form Ini = [x10,dx10dt x20 dx20dt x30 dx30dt x40 dx40dt x50 dx50dt x60 dx60dt]

close all
Ini = [0,0,0,0,0,0,0.1,0,0,0,0,0];

theta0 = 0;

dxdt = zeros(6,n_int);
dxdt(:,1) = Ini(2:2:end);
x = zeros(6,n_int);
x(:,1) = Ini(1:2:end);
dtheta = 0;
theta = zeros(n_int,1);
theta = theta0;

for i=2:n_int
    t_int(i-1) = (i-2)*deltaT;
    inp = num2cell(x(:,i-1));
    dxdt(:,i) = dxdt(:,i-1) + acc(theta(i-1),inp{:})*deltaT;
    x(:,i) = x(:,i-1) + dxdt(:,i)*deltaT;
    theta(i) = theta(i-1) + dtheta*deltaT;
end
t_int(i) = (i-1)*deltaT;

leg = {};

X = sym(zeros(6,1));
X(:) = [x1,x2,x3,x4,x5,x6];

for i=1:size(x,1)
    if i==6
        style = '--';
    else
        style = '-';
    end
    plot(t_int,x(i,1:length(t_int)),'linewidth',3,'linestyle',style)
    hold on
    leg{i} = string(X(i));
end
legend(leg{:},'interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('Displacement [m]','interpreter','latex');
set(gca,'fontsize',28);

input('Press enter to animate');
animateFlexibleStructure(x,theta,t_int)








