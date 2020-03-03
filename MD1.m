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
ITA1 = [-4*BeamStiff_beam(L1,I1)*x1(t), 0, 0].';
IRA1 = [0,R1,0].';

ITA2 = [4*BeamStiff_beam(L2,I1)*(x2(t)-x1(t)),0,0].';
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
kb = BeamStiff_blade(L5,I2);


K = [2*k1+2*k2  -2*k2      0         0         0   0;
     -2*k2 2*k2+2*k3  -2*k3     0         0   0;
     0     -2*k3      2*k3+2*k4 -2*k4     0   0;
     0     0          -2*k4     2*k4+2*kb -kb -kb;
     0     0          0         -kb       kb  0;
     0     0          0         -kb       0   kb];
nulmat = zeros(size(M));

A = [M,nulmat;
     nulmat,M];
B = [nulmat,K;
     -M,nulmat]; 
 
[u,lambda] = eig(-B,A);

lambda = diag(lambda);

[lambda,ind] = sort(lambda);
omega0 = lambda./(2*pi);
u = u(:,ind);
phi = u(1:6,:);
PlotModeShapes(phi,omega0);

keyboard


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

[A,b] = equationsToMatrix(eqsys(index),dx.');
syms xx1 xx2 xx3 xx4 xx5 xx6

b = subs(b,{x1,x2,x3,x4,x5,x6},{xx1,xx2,xx3,xx4,xx5,xx6});
[C,d] = equationsToMatrix(b,[xx1,xx2,xx3,xx4,xx5,xx6]);

keyboard

% numerical solution
acc = linsolve(A,b);
acc = acc(1:6);


u = sym(zeros(6,1));
u = [u1 u2 u3 u4 u5 u6];
acc = subs(acc,[x1,x2,x3,x4,x5,x6,theta(t)],[u(1),u(2),u(3),u(4),u(5),u(6),phi]);


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








