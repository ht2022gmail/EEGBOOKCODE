function draw_cov_ellipse(S, opt)

[V, D] = eig(S);
theta = linspace(0,2*pi,100);

x0 = 2*[D(1,1)*cos(theta); D(2,2)*sin(theta)];
x = V*x0;
if nargin<2
    plot(x(1,:), x(2,:), 'k--', 'linewidth', 2);
else
    plot(x(1,:), x(2,:), opt, 'linewidth', 2);
end
 
