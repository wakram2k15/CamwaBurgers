function [L2e,H1e]=Err_gaussiantest(coord1,coord,elem,ed2el,n2ed,temp1,uh_1)
L2e_1=0; H1e_1=0; 

nt=size(elem,1);      % size of number of element in the previous refinement
ne=size(ed2el,1);     % size of number of edges in the previous refinement
marker=sparse(ne,1);   
marker(:,1)=[(size(coord,1)+1):size(coord1,1)]';
for i=1:nt,           %loop over the elements in previous refinement 
    ct=elem(i,:);     %nodes of current elements in previous refinement
    curcoords=coord(ct,:);   %coordinates of nodes of current element in previous refinement
    Cu1=temp1(ct);        % Q11 approximate solution at nodes of current elements
    
    %=================constructing the 4 triangles in current refinement from 
    %====================current element in previous refinement============
    ce=diag(n2ed(ct([2 3 1 ]),ct([3 1 2])));
    m1=marker(ce(1));  
    m2=marker(ce(2));
    m3=marker(ce(3));
    curelem=zeros(4,3);
    curelem(1,:)=[m1 m2 m3];
    curelem(2,:)=[ct(1) m3 m2]; 
    curelem(3,:)=[ct(2) m1 m3]; 
    curelem(4,:)=[ct(3) m2 m1];
for j=1:4                                      %loop over those 4 triangles
       cur_nodes=curelem(j,:);                 % nodes of that triangle in current refinement
       cord_curelem =coord1(cur_nodes,:);%coordinates of nodes of that triangle in current refinement
       p1=cord_curelem(1,:); p2=cord_curelem(2,:); p3=cord_curelem(3,:);
    
       x1 = p1(1); y1 = p1(2);
       x2 = p2(1); y2 = p2(2);
       x3 = p3(1); y3 = p3(2);
       xw = TriGaussPoints(5);                % Gauss quadrature points             
       Area = abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;
       NP=length(xw(:,1));
       z1 = 0; 
                                  % evaluated at Gauss quadrature points
for k = 1:NP                      % loop over quadrature points
     x = x1*(1-xw(k,1)-xw(k,2))+x2*xw(k,1)+x3*xw(k,2);
     y = y1*(1-xw(k,1)-xw(k,2))+y2*xw(k,1)+y3*xw(k,2);

    cord = [x y];
 z1 = z1 + (temp_curcg(cord,Cu1,curcoords)-temp_curcg(cord,uh_1(cur_nodes),cord_curelem))^2*xw(k,3);
end
  
  L2e_1 = L2e_1 + Area*z1;         %L2e_1 is L2 error corresponding to first component Q11
end
    
    
   P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mp12=1/2*(P1+P2); mp13=1/2*(P1+P3); mp23=1/2*(P3+P2);

    mp12_nodes=find(coord1(:,1)==mp12(1,1) & coord1(:,2)==mp12(1,2));
    mp13_nodes=find(coord1(:,1)==mp13(1,1) & coord1(:,2)==mp13(1,2));
    mp23_nodes=find(coord1(:,1)==mp23(1,1) & coord1(:,2)==mp23(1,2));

    mk=1/2*det([1 P1;1 P2;1 P3]);   % mk is area of triangle
      L1=[1 mp12;1 mp23;1 mp13]\[1 0 0]'; %barycentric coordinates corresponding to triangle in current refinement
      L2=[1 mp12;1 mp23;1 mp13]\[0 1 0]';
      L3=[1 mp12;1 mp23;1 mp13]\[0 0 1]';
      
    uxh1=uh_1(mp12_nodes)*L1(2)+uh_1(mp23_nodes)*L2(2)+uh_1(mp13_nodes)*L3(2);
    uyh1=uh_1(mp12_nodes)*L1(3)+uh_1(mp23_nodes)*L2(3)+uh_1(mp13_nodes)*L3(3);
    
    
    % barycentric coordinates corresponding to triangle in previous refinement
      L1=[1 P1;1 P2;1 P3]\[1 0 0]';
      L2=[1 P1;1 P2;1 P3]\[0 1 0]';
      L3=[1 P1;1 P2;1 P3]\[0 0 1]';
      
    temp1xh1=Cu1(1)*L1(2)+Cu1(2)*L2(2)+Cu1(3)*L3(2);
    temp1yh1=Cu1(1)*L1(3)+Cu1(2)*L2(3)+Cu1(3)*L3(3);
    
    H1e_1=H1e_1+mk*(temp1xh1-uxh1)^2+mk*(temp1yh1-uyh1)^2;
   
end
 L2e=sqrt(L2e_1);%/Norm(coord, elem, temp1); 
 H1e=sqrt(H1e_1);%/H1_Norm(coord, elem, temp1);
 
end


