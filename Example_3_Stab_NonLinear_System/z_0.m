function uv=z_0(z)
x=z(1); y=z(2); t=z(3);
% uv=sin(pi*x)*sin(pi*y)- x*y*(1 - x)*(1 - y);
uv=x*y*(1 - x)*(1 - y)- sin(pi*x)*sin(pi*y);