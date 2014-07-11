theta = linspace(0,90,100);
r = linspace(1,100,100);
r_0 = r+1;
v = linspace(1,10,100);
omega = linspace(1,100,100);
A = 50;
R = 5;
b=1
a=1
g00 = 1;
delta = 1;
gij = [1,1,1]
n_a = 3.0;
n_c = 0.8;
n_d = 0.2;
sigma = 1;

%refractive index for a CIPT
fun = @(r) exp(-r.^2).*log(r).^2./r;
chi = @(r) -2.*r.^2;
q = integral(fun,1,inf);
nCIPT = A./r.*exp(q);

% calculate the boradband refractive index for circularly sym photonic
% blackhole here
h_i = sqrt(gij)
h1 = sqrt(1)
h2 = sqrt(1)
h3 = sqrt(1)

g = g00*((b./r).^2+(1-a./r).^2)
eps = delta*(h1*h2*h3)./(h_i*sqrt(g00))
mu = eps
nPBH = sqrt(g/g00)

n = nCIPT
%potential or kinetic energy in the total zero energy state 
U = -n.^2/2

% gaussian distribution
n1 = n_a.*exp(-r.^2./sigma.^2)+n_c;

%mexican hat distribution
n2 = (n_a-(n_d.*r.^2)/(sigma.^2)).*exp(-r.^2./sigma.^2)+n_c;

% double attractors
n3 = n_a*(exp(-(r-r_0).^2./sigma.^2)+exp(-(r-r_0).^2./sigma.^2)) + n_c;

% g = [1,1;0,1];
% dx_dtau = [1,1;1,1];
% dy_dtau = [1,1;1,1];
% 
% %see  alegri et al eq1
% % lagrangian
% L = (1/2)*g.*dx_dtau.*dy_dtau

%eq 7-10
% ray tracing 

% r_dd = r.*phi_d.^2+t_d.^2*delta_r*n./n.^3;
% phi_dd = (-2*r*r_d*phi_d+(t_d.^2*d_phi*n)/n.^3))./r.^2;
% phi_dd = t_d*delta_z*n./n.^3;
% t_dd = (t_d*(t_d*delta_t*n+2*(z_dd*delta_z*n+phi_d*delta_phi*n+phi_d*delta_phi*n+r_d*delta_r*n))./n;

%reduction of ray tracing equations
%omega = diff(theta);
%a = diff(diff(r);
%v = diff(r);
deltar_n = n;
deltatheta_n = linspace(1,1,100);
deltaomegan = linspace(1,1,100);

%a = n.*r.*omega.^2-2.*v.*omega.*deltatheta_n-v.^2.*deltar_n+r.^2.*omega.^2.*deltar_n./n;
%alpha = -2.*n.*r.*v.*omega-v.^2.*deltaomegan+r.^2.*omega.^2.*deltaomegan+2.*r.^2.*v.*omega.deltar_n./(n.*r.^2);

%plot

[t,r] = meshgrid(theta,nCIPT);
[x,y] = pol2cart(t,r);