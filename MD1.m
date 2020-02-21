clear all
syms E I1 I2 L L1 L2 L3 L4 L5 L6 m1 m2 m3 m4 m5 m6 x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) phi(x) t g theta(t) R1(t) R2(t) R3(t) R4(t) R5(t) R6(t)

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
b3Ac = b2Ab + diff(b2rbc,t);
b4Ad = b3Ac + diff(b3rcd,t);

%Accelerations for blade tips:
b5Ae = T_theta*b4Ad + cross(diff(omega,t),b5rde) + cross(omega,cross(omega,b5rde)) + 2 * cross(omega,diff(b5rde,t)) + diff(b5rde,t,t);
b5Af = T_theta*b4Ad + cross(diff(omega,t),b5rdf) + cross(omega,cross(omega,b5rdf)) + 2 * cross(omega,diff(b5rdf,t)) + diff(b5rdf,t,t);



BeamStiff = @(L,I) 3*E*I/(L^3)

IPA = [0, -m1*g, 0].';
ITA1 = [-BeamStiff(L1,I1)*x1(t), 0, 0].';
IRA1 = [0,R1,0].';

ITA2 = [BeamStiff(L2,I1)*(x2(t)-x1(t))].';
IRA2 = [0,-R2,0].';

IPB = [0, -m2*g, 0].';
ITB2 = -ITA2;
IRB2 = -IRA2;

ITB3 = [(x3(t)-x2(t))*BeamStiff(L3,I1),0,0].';
IRB3 = [0,-R3,0].';

IPC = [0,-m3*g,0].';
ITC3 = -ITB3;
IRC3 = -IRB3;
ITC4 = [(x4(t)-x3(t))*BeamStiff(L4,I1),0,0].';
IRC4 = [0,-R4,0].';

IPD = [0,R4,0].';
ITD4 = -ITC4;
IRD4 = -IRC4;
%

IRD5 = [0,-R5,0].';

IPE = T_theta*[0,-m5*g,0].';
ITE = T_theta*[-BeamStiff(L5,I2)*x5(t),0,0].';
IRE = [0,R5,0].';



IPF = T_theta*[0,-m6*g,0].';
ITF = T_theta*[-BeamStiff(L6,I2)*x6(t),0,0].';
IRF = T_theta*[0,R6,0].';

ITD5 = T_theta.'*(ITE+ITF)

eq1(:) = m1.*(b1Aa) - ( IPA+IRA1+IRA2+ITA1+ITA2);
eq2(:) = m2.*(b2Ab) - ( IPB + IRB2 + IRB3 + ITB2 + ITB3);
eq3(:) = m3.*(b3Ac) - ( IPC + IRC3 + IRC4 + ITC3 + ITC4);
eq4(:) = m4.*(b4Ad) - ( IPD + IRD4 + IRD5 + ITD4);
eq5(:) = m5.*(b5Ae) - ( IPE + IRE + ITE - ITD5);
eq6(:) = m6.*(b5Af) - ( IPF + IRF + ITF - ITD5);

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

b2 = 0.0350;
h2 = 0.0015;
I2 = b2*h2^3/12;

g = 9.81;

eq1 = subs(eq1);
eq2 = subs(eq2);
eq3 = subs(eq3);
eq4 = subs(eq4);
eq5 = subs(eq5);
eq6 = subs(eq6);

eqsys = [eq1;eq2;eq3;eq4;eq5;eq6];
syms ddx1 ddx2 ddx3 ddx4 ddx5 ddx6 RA RB RC RD RE RF
dderivatives = [ddx1 ddx2 ddx3 ddx4 ddx5 ddx6];
forces = [RA RB RC RD RE RF]
eqsys = subs(eqsys,[diff(x1,t,t),diff(x2,t,t),diff(x3,t,t),diff(x4,t,t),diff(x5,t,t),diff(x6,t,t),R1,R2,R3,R4,R5,R6],[dderivatives,forces])
dx = [dderivatives,forces];
dderivatives = [ddx1 ddx2 ddx3 ddx4 ddx5 ddx6];


index = [];

for i=1:size(eqsys,1)
    if mod(i,3) ==0
        
    else
    index(end+1) = i;
    end
end

keyboard
[A,b] = equationsToMatrix(eqsys(index),dx)


[V,S] = odeToVectorField([eq1;eq2;eq3;eq4;eq5]);
M = matlabFunction(V,'vars',{'t','Y'});

%Define mode shapes
Beta1 = 1.87510407/L;
Sigma1 = 0.7341;
Phi1(x) = cosh(Beta1*(L-x))+cos(Beta1*(L-x))-Sigma1*(sinh(Beta1*(L-x))+sin(Beta1*(L-x)));

y1(t) = q1(t)*Phi1(L);
y2(t) = q2(t)*Phi1(L);
y3(t) = q3(t)*Phi1(L);
y4(t) = q4(t)*Phi1(L);
y5(t) = q5(t)*Phi1(L);
y6(t) = q6(t)*Phi1(L);


eq1 = subs(eq1);
eq2 = subs(eq2);
eq3 = subs(eq3);
eq4 = subs(eq4);
eq5 = subs(eq5);
eq6 = subs(eq6);

eqsys = [eq1;eq2;eq3;eq4;eq5;eq6];

eqsys = vpa(eqsys);



[V,S] = odeToVectorField([eq1;eq2;eq3;eq4;eq5;]);
M = matlabFunction(V,'vars',{'t','Y'});

%Define initial conditions
Ini = [0,0.1,0,0,0,0,0,0,1,0];

opts = odeset('Stats','on');

simTime = 15;
disp('Solving');

%Solve
sols = ode45(M,[0,simTime],Ini,opts);

for i=1:length(S)
    plot(sols.x,sols.y(i,:),'linewidth',3)
    leg{i} = string(S(i));
    hold on
end
legend(leg{:});
















