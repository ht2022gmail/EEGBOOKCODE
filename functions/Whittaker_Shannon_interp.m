function y = Whittaker_Shannon_interp(t, dt, t0, x0)

A = repmat(t', 1, length(t0)) - ones(length(t), 1)*t0;
A = sinc(A/dt);

y = x0*A';