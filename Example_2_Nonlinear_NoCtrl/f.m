function fv = f(z)
x = z(1); y = z(2);  t = z(3);
shift=0; eta=5;


% %%z=exp(-t)sin(pix)sin(piy)
% fv=(-1+2*pi^2*eta)*exp(-t)*sin(pi*x)*sin(pi*y) + (x*(1 - x)* y *(1 - y) + exp(-t)*sin(pi*x)*sin(pi*y))*...
%     exp(-t)*pi*(cos(pi*x)*sin(pi*y)+sin(pi*x)*cos(pi*y))+...
%     ((1-2*x)*y*(1-y)+x*(1-x)*(1-2*y))*exp(-t)*sin(pi*x)*sin(pi*y);



%%z=exp(t)sin(pix)sin(piy)
fv=(1+2*pi^2*eta)*exp(t)*sin(pi*x)*sin(pi*y) + (x*(1 - x)* y *(1 - y) + exp(t)*sin(pi*x)*sin(pi*y))*...
    exp(t)*pi*(cos(pi*x)*sin(pi*y)+sin(pi*x)*cos(pi*y))+...
    ((1-2*x)*y*(1-y)+x*(1-x)*(1-2*y))*exp(t)*sin(pi*x)*sin(pi*y);