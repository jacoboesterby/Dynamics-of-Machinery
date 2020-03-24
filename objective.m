function val = objective(x)
% x = [m1,m2,m3,m4,m5,m6,k1,k2,k3,k4,k5,k6];

xi = [0.0024,0.0018,0.0020,0.0031,0.0028,0.0040].';

M = [x(1),0,0,0,0,0;
    0,x(2),0,0,0,0;
    0,0,x(3),0,0,0;
    0,0,0,x(4),0,0;
    0,0,0,0,x(5),0;
    0,0,0,0,0,x(6)];

K = [2*x(7)+2*x(8)  -2*x(8)      0         0         0   0;
    -2*x(8) 2*x(8)+2*x(9)  -2*x(9)     0         0   0;
    0     -2*x(9)      2*x(9)+2*x(10) -2*x(10)     0   0;
    0     0          -2*x(10)     2*x(10)+x(11)+x(12) -x(11) -x(12);
    0     0          0         -x(11)       x(11)  0;
    0     0          0         -x(12)       0   x(12)];
nulmat = zeros(size(M));

alpha = 0; % 2.296811617444550e-02;
beta = 0; % 1.675208490914375e-04;

D3 = alpha.*M+beta.*K;

A = [M,nulmat;
    nulmat,M];
B = [nulmat,K;
    -M,nulmat];

[u,lambda] = eig(-B,A);

omegaex = [0.84;1.22;1.77;3.64;5.45;7.26];

lambda = diag(lambda);
[lambda,ind] = sort(lambda);
omega0 = imag(lambda);

for k=1:6
    lambda((2*k-1):2*k) = [-xi(k)*abs(omega0(2*k))-abs(omega0(2*k))*sqrt(1-xi(k)^2)*i;-xi(k)*abs(omega0(2*k))+abs(omega0(2*k))*sqrt(1-xi(k)^2)*i];
    omegad(k) = omega0(2*k)*sqrt(1-xi(k)^2)/(2*pi);
end
val = norm(abs(omegad.' - omegaex));
end