function [ysxv,ysyv]=ysxe(z)
x=z(1); y=z(2); t=z(3);
% ysxv=(1 - 2*x)*y*(1 - y);
% ysyv=(1 - 2*y)*x*(1 - x);

ysxv=pi * cos(pi*x) * sin(pi*y);
ysyv=pi*sin(pi*x) * cos(pi*y);