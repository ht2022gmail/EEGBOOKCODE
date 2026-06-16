function [y, W] = infomaxICA(x)
[n, m] = size(x);

eta = 0.01;
Wold = eye(n);
Niter = 10000;
err = 1;

iter = 1;
while iter<Niter && err>10e-6
    y = Wold*x;
    phi = tanh(y);
    Wnew = Wold + eta*(eye(n)-phi*y'/m)*Wold;
    err = norm(Wnew-Wold,'fro')/norm(Wold,'fro');
    Wold = Wnew;
    iter = iter+1;
end

y = Wnew*x;
W = Wnew;