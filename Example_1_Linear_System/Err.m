function [L2e,H1e]=Err(Coord, Elem, uh, u,T)

L2e=0; H1e=0; 

for j=1:size(Elem,1),
    curnodes=Elem(j,:);
    curcoords=Coord(curnodes,:);
    Cuh=uh([curnodes,1]);
    Cu=u([curnodes,1]);
    P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mp12=1/2*(P1+P2); mp13=1/2*(P1+P3); mp23=1/2*(P3+P2);
    cg=1/3*(P1+P2+P3);
    
    uhmp12=1/2*(Cuh(1)+Cuh(2)); uhmp13=1/2*(Cuh(1)+Cuh(3));
    uhmp23=1/2*(Cuh(3)+Cuh(2)); uhcg=1/3*(Cuh(1)+Cuh(2)+Cuh(3));
    
    mk=1/2*det([1 P1;1 P2;1 P3]);
    
    L2e=L2e+mk/60*(8*(ue([mp12,T])-uhmp12)^2+8*(ue([mp13,T])-uhmp13)^2+8*(ue([mp23,T])-uhmp23)^2+...
            3*(Cu(1)-Cuh(1))^2+3*(Cu(2)-Cuh(2))^2+3*(Cu(3)-Cuh(3))^2+27*(ue([cg,T])-uhcg)^2);

    L1=[ones(1,3);curcoords']'\[1;0;0];
    L2=[ones(1,3);curcoords']'\[0;1;0];
    L3=[ones(1,3);curcoords']'\[0;0;1];
    uxh=Cuh(1)*L1(2)+Cuh(2)*L2(2)+Cuh(3)*L3(2);
    uyh=Cuh(1)*L1(3)+Cuh(2)*L2(3)+Cuh(3)*L3(3);
   
    [uxev uyev]=uxe([cg,T]);
    
    H1e=H1e+mk*(uxev-uxh)^2+mk*(uyev-uyh)^2;
   end
   L2e=sqrt(L2e); H1e=sqrt(H1e);
end







