r = linspace(1,100,100);
A = 50;
R = 5;
fun = @(r) exp(-r.^2).*log(r).^2./r;
q = integral(fun,1,inf);
n = A./r.*exp(q);
