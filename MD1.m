clear all
syms L L1 L2 L3 L4 L5 L6 m1 m2 m3 m4 m5 m6 y1(t) y2(t) y3(t) y4(t) y5(t) y6(t) x1 x2 x3 x4 x5 x6 phi(x) t g theta(t) q1(t) q2(t) q3(t) q4(t) q5(t) q6(t)

%Position vectors
iroa = [y1(t),L1,0].';
b1rab = [y2(t),L2,0].';
b2rbc = [y3(t),L3,0].';
b3rcd = [y4(t),L4,0].';
b5rde = [y5(t),L5,0].';
b5rdf = [y6(t),-L6,0].';

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

keyboard
%Define gravitational force:
g = [0,-g,0].';

eq1 = m1.*(b1Aa+g);
eq2 = m2.*(b2Ab+g);
eq3 = m3.*(b3Ac+g);
eq4 = m4.*(b4Ad+g);
eq5 = m5.*(b5Ae+T_theta*g);
eq6 = m6.*(b5Af+T_theta*g);


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

Ebeam = 2*10^11;
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

index = find(eqsys~=0);

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
















