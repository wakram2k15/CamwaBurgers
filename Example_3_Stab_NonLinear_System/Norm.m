function L2e=Norm(Coord, Elem, uh)

L2e=0; 

for j=1:size(Elem,1),
    curnodes=Elem(j,:);
    curcoords=Coord(curnodes,:);
    Cuh=uh([curnodes,1]);
%     Cu=u([curnodes,1]);
    P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mp12=1/2*(P1+P2); mp13=1/2*(P1+P3); mp23=1/2*(P3+P2);
    cg=1/3*(P1+P2+P3);
    
    uhmp12=1/2*(Cuh(1)+Cuh(2)); uhmp13=1/2*(Cuh(1)+Cuh(3));
    uhmp23=1/2*(Cuh(3)+Cuh(2)); uhcg=1/3*(Cuh(1)+Cuh(2)+Cuh(3));
    
    mk=1/2*det([1 P1;1 P2;1 P3]);
    
    L2e=L2e+mk/60*(8*(uhmp12)^2+8*(uhmp13)^2+8*(uhmp23)^2+...
            3*(Cuh(1))^2+3*(Cuh(2))^2+3*(Cuh(3))^2+27*(uhcg)^2);

    
   end
   L2e=sqrt(L2e); 
end







