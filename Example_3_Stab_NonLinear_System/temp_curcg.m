function T=temp_curcg(z,temp,curcoord)
x=z(1);
y=z(2);
Q1=curcoord(1,:); Q2=curcoord(2,:); Q3=curcoord(3,:);
L1=[1 Q1;1 Q2;1 Q3]\[1 0 0]';
L2=[1 Q1;1 Q2;1 Q3]\[0 1 0]';
L3=[1 Q1;1 Q2;1 Q3]\[0 0 1]';
lamda1=L1(1)+L1(2)*x+L1(3)*y;
lamda2=L2(1)+L2(2)*x+L2(3)*y;
lamda3=L3(1)+L3(2)*x+L3(3)*y;
T=temp(1)*lamda1 + temp(2)*lamda2 +temp(3)*lamda3;
end