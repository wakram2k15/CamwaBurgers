function [ST,Mass,Ah1,Nh,V]=stima_all(vertices,z)
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

%Calculate matrix with mix term around non-zero velocity
Ah1=(1/3)*mk*[L1(2)+L1(3) L1(2)+L1(3) L1(2)+L1(3)
        L2(2)+L2(3) L2(2)+L2(3) L2(2)+L2(3)
        L3(2)+L3(3) L3(2)+L3(3) L3(2)+L3(3)];

Nh=(z(1)*(L1(2)+L1(3))+z(2)*(L2(2)+L2(3))+z(3)*(L3(2)+L3(3)))*Mass;
%Mass=(1/3)*mk*[1 0 0;0 1 0;0 0 1];

%%%caculate Jacobian matrix
a1=L1(2)+L1(3); a2=L2(2)+L2(3); a3=L3(2)+L3(3);
b1=z(1)/6+z(2)/12+z(3)/12; b2=z(1)/12+z(2)/6+z(3)/12; b3=z(1)/12+z(2)/12+z(3)/6;

V=[a1*b1 a1*b2 a1*b3 ; a2*b1 a2*b2 a2*b3 ; a3*b1 a3*b2 a3*b3];
