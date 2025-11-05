clear all;
format long;

% Increase nl for more refine case

[Coord,Elem,Nb,Db]=InitialMesh(1);


T=0.1; shift=24; 
eta=1;

for nl=1:6
   
%% Edge-Node-Element Connections
[n2ed,ed2el]=edge(Elem,Coord);

%% Element Redrefine
[Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);

%%%%%%%%%%% find h %%%%%%%%%%%%%%
h(nl)=sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
dt=0.001;

FullNodes=[1:size(Coord,1)]; % Global node numbers
FreeNodes=setdiff(FullNodes, unique(Db));
    nfull=length(FullNodes);
    nfree=length(FreeNodes);

% Intializing the matrices
M_1=sparse(size(Coord,1),size(Coord,1)); 
K_1=sparse(size(Coord,1),size(Coord,1));
Ah2=sparse(size(Coord,1),size(Coord,1));
Ah3=sparse(size(Coord,1),size(Coord,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly of A 
% stima is element stiffness matrices
%[]
    for j=1:size(Elem,1)
        [ST,Mass,Mix]=stima(Coord(Elem(j,:),:));
        K_1(Elem(j,:),Elem(j,:))=K_1(Elem(j,:),Elem(j,:))+ST;
        M_1(Elem(j,:),Elem(j,:))=M_1(Elem(j,:),Elem(j,:))+Mass; 
        Ah2(Elem(j,:),Elem(j,:))=Ah2(Elem(j,:),Elem(j,:))+...
            ys([sum(Coord(Elem(j,:),:))/3,0])*Mix; 
        [ysxev ysyev]=ysxe([sum(Coord(Elem(j,:),:))/3,0]);
        Ah3(Elem(j,:),Elem(j,:))=Ah3(Elem(j,:),Elem(j,:))+...
            (ysxev+ysyev)*Mass;
    end
    
    A_h=-eta*K_1(FreeNodes,FreeNodes)+shift*M_1(FreeNodes,FreeNodes)-...
        Ah2(FreeNodes,FreeNodes)-Ah3(FreeNodes,FreeNodes);
    M_h=M_1(FreeNodes,FreeNodes);
    
    yh=zeros(length(FullNodes),1);
    
    y_temp=zeros(length(FullNodes),1);
        
    for j=1:size(Elem,1),
                y_temp(Elem(j,:))=y_temp(Elem(j,:))+...
                    det([1 1 1; Coord(Elem(j,:),:)'])*y_0([sum(Coord(Elem(j,:),:))/3,0])/6-...
                 det([1 1 1; Coord(Elem(j,:),:)'])*ys([sum(Coord(Elem(j,:),:))/3,0])/6;
    end 

        
    y_temp=M_h\y_temp(FreeNodes); 
    

   uh=zeros(length(FullNodes),1);
   for k=1 %BDF1 in first time step
       y_temp1=(M_h-dt*A_h)\((M_h)*y_temp);
       yh(FreeNodes)=y_temp1;
       E(k,nl)=Norm(Coord, Elem, yh);
       Eu(k,nl)=Norm(Coord, Elem, uh);
   end
%     
   for k=2:T/dt %BDF2 in second time step onwards
       yh(FreeNodes)=(1.5*M_h-dt*(A_h))\(M_h*(2*y_temp1-0.5*y_temp));
       y_temp=y_temp1;
       y_temp1=yh(FreeNodes);
       %yh(FreeNodes)=yh(1:nfree,1);
       E(k,nl)=Norm(Coord, Elem, yh);%+Norm(Coord, Elem, zh);    
       Eu(k,nl)=Norm(Coord, Elem, uh);
   end

    if (nl > 1)
        [L2e(nl-1),H1e(nl-1)]= ...
        Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempy, yh);
        [L2eu(nl-1),H1eu(nl-1)]= ...
        Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempu, uh);
    end
        
        tempy=yh;
        tempu=uh;
        Coord1=Coord;
        Elem1=Elem;
    
 
end


for ml=1:(nl-2)
    ocl2(ml+1)=log(L2e(ml)/L2e(ml+1))/log(h(ml)/h(ml+1));
    och1(ml+1)=log(H1e(ml)/H1e(ml+1))/log(h(ml)/h(ml+1));
end
ocl2;
och1;

oczl2(1)  = NaN;
oczh1(1)  = NaN;
format long e
Order_of_Convergence=table(L2e', ocl2', H1e', och1');
Order_of_Convergence.Properties.VariableNames = {'L2_Err_z' 'L2_Conv_z'...
                'H1_Err_z' 'H1_Conv_z'}


%%%Log-log plot
figure;
loglog(h(2:end), L2e, 's-', 'LineWidth', 1.1, 'MarkerSize', 8);
hold on;
loglog(h(2:end), H1e, 's-', 'LineWidth', 1.1, 'MarkerSize', 8);
grid on;

xlabel('Discretization parameter h');
ylabel('Errors');
title('Log–Logplot of errors against discretization parameter h');
legend('err_{L^2}(z_{h_i})', 'err_{H^1}(z_{h_i})');

