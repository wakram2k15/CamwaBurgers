function uv=ye(z)
x=z(1); y=z(2); t=z(3);
uv = sin(pi*x)*sin(pi*y)-x*y*(x - 1)*(y - 1);
% uv = exp(-2*pi^2*t/5)*sin(pi*x)*sin(pi*y);
%uv=exp(-t)*sin(pi*x)*sin(pi*y);
%uv=exp(40*t)*sin(pi*x)*sin(pi*y);
%uv=exp((40-2*pi^2)*t)*sin(pi*x)*sin(pi*y);
%uv=exp(t)*sin(pi*x)*sin(pi*y);