function [ST,Mass]=stima(vertices)
% area of current element 
mk=1/2*det([ones(1,3);vertices']);

% barycentic coordinates of current element
L1=[ones(1,3);vertices']'\[1;0;0];
L2=[ones(1,3);vertices']'\[0;1;0];
L3=[ones(1,3);vertices']'\[0;0;1];

% element stiffness matrix
ST=mk*[L1(2)*L1(2)+L1(3)*L1(3) L1(2)*L2(2)+L1(3)*L2(3) L1(2)*L3(2)+L1(3)*L3(3)
       L2(2)*L1(2)+L2(3)*L1(3) L2(2)*L2(2)+L2(3)*L2(3) L2(2)*L3(2)+L2(3)*L3(3)
       L3(2)*L1(2)+L3(3)*L1(3) L3(2)*L2(2)+L3(3)*L2(3) L3(2)*L3(2)+L3(3)*L3(3)];
   
Mass=mk*[1/6 1/12 1/12;1/12 1/6 1/12;1/12 1/12 1/6]; % mass matrix
%Mass=(1/3)*mk*[1 0 0;0 1 0;0 0 1];
