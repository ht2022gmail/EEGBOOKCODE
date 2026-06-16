function dydt = JansenRit_singlecolumnmodel(t, y, pt, p)

p = interp1(pt, p, t); % Interpolate the data set (pt,p) at time t

% A, B, a, C1, C2, C3, C4
A = 3.25;
B = 22;
C = 135;
a = 100;
b = 50;
C1 = C;
C2 = 0.8*C;
C3 = 0.25*C;
C4 = 0.25*C;


dydt = zeros(6,1);
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);
dydt(4) = A*a*sigmoid1(y(2)-y(3)) - 2*a*y(4) - a^2*y(1);
dydt(5) = A*a*(p + C2*sigmoid1(C1*y(1))) - 2*a*y(5) - a^2*y(2);
dydt(6) = B*b*(0.3*p + C4*sigmoid1(C3*y(1))) - 2*b*y(6) - b^2*y(3);


function y = sigmoid1(v)
% e0, r, v0
e0 = 2.5;
r = 0.56;
v0 = 6;

y = 2*e0./(1+exp(r*(v0-v)));