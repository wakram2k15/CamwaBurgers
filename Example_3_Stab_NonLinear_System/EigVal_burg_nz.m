clear all;
format long;

% Geometry of finite elements - triangulation, local global node numbering,
% coordinates, boundary nodes

[Coord,Elem,Nb,Db]=InitialMesh(1);


shift=3; 
%eta=1;
eta =0.1;

for nl=1:3
   
%% Edge-Node-Element Connections
[n2ed,ed2el]=edge(Elem,Coord);

%% Element Redrefine
[Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);

%%%%%%%%%%% find h %%%%%%%%%%%%%%
h(nl)=sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);


FullNodes=[1:size(Coord,1)]; % Global node numbers
FreeNodes=setdiff(FullNodes, unique(Db));
    nfull=length(FullNodes);
    nfree=length(FreeNodes);

% Intializing the matrices
M=sparse(size(Coord,1),size(Coord,1)); 
K=sparse(size(Coord,1),size(Coord,1));
A1=sparse(size(Coord,1),size(Coord,1));
A2=sparse(size(Coord,1),size(Coord,1));



% Assembly of A 
% stima is element stiffness matrices
%[]
    for j=1:size(Elem,1)
        [ST,Mass,Ah1]=stima_1(Coord(Elem(j,:),:));
        K(Elem(j,:),Elem(j,:))=K(Elem(j,:),Elem(j,:))+ST;
        M(Elem(j,:),Elem(j,:))=M(Elem(j,:),Elem(j,:))+Mass; 
        A1(Elem(j,:),Elem(j,:))=A1(Elem(j,:),Elem(j,:))+...
            ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
        [ysxev ysyev]=ysxe([sum(Coord(Elem(j,:),:))/3,0]);
        A2(Elem(j,:),Elem(j,:))=A2(Elem(j,:),Elem(j,:))+...
            (ysxev+ysyev)*Mass;
    end
    
    Ah=-eta*K(FreeNodes,FreeNodes)+shift*M(FreeNodes,FreeNodes)-...
        A1(FreeNodes,FreeNodes)-A2(FreeNodes,FreeNodes);
    Mh=M(FreeNodes,FreeNodes);
    
    t=length(FreeNodes);
    [V,D]=eigs(Ah,Mh,t);

    j=1;
    for i=1:t-1
        if real(D(i,i))>0% && real(D(i))<4
            UE(nl,j)=D(i,i);
            j=j+1;
        end
    end

 
end


plot(real(D),imag(D),'b*');
%axis([-100 30 -10 10])

for k=1:nl-1
    Err1(k)=abs(UE(k,1)-UE(k+1,1));
end

for ml=1:(nl-2)
    ocEig1(ml)=log(Err1(ml)/Err1(ml+1))/log(h(ml)/h(ml+1));
end
ocEig1









