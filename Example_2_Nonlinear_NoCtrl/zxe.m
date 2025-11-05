function [uxv,uyv]=zxe(z)
x=z(1); y=z(2); t=z(3);

% uxv=pi*exp(-t)*cos(pi*x)*sin(pi*y);
% uyv=pi*exp(-t)*sin(pi*x)*cos(pi*y);

uxv=pi*exp(t)*cos(pi*x)*sin(pi*y);
uyv=pi*exp(t)*sin(pi*x)*cos(pi*y);