function [y, yb, tb] = blocky2(X, w, t, tau)
% computation of 1st TRC (y)
y = w(:,1)'*X;   
y = y-mean(y); y = y/std(y);
% y = sign(sum(y*X'))*y;  % this ensures the correct sign.

% blocking of y
tb = repmat(t', 1, tau) + repmat(0:tau-1, length(t),1);
yb = y(tb);
yb = bsxfun(@minus, yb, mean(yb,2));

