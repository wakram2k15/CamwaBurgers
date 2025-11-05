clear all;
format long;

[Coord,Elem,Nb,Db]=InitialMesh(1);

T=0.1;
shift=0; eta=5;

for nl=1:5
   
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


y_temp=zeros(length(FullNodes),1);
    for j=1:length(FullNodes)
       y_temp(FullNodes(j),1) = y_0([Coord(FullNodes(j),:),0]);
    end
y_temp=y_temp(FreeNodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly of A 

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


        for k=1 %BDF1 in first time step
            b=zeros(size(Coord,1),1); % global load vector
            for j=1:size(Elem,1),
                b(Elem(j,:))=b(Elem(j,:))+det([1 1 1; Coord(Elem(j,:),:)'])*...
                    f([sum(Coord(Elem(j,:),:))/3,k*dt])/6;
            end   
              b=b(FreeNodes);
        
            y_temp1=(M_h-dt*A_h)\(M_h*y_temp+dt*b);
            yh(FreeNodes)=y_temp1;
        end
%     
        for k=2:T/dt %BDF2 in second time step onwards
             b=zeros(size(Coord,1),1); % global load vector
            for j=1:size(Elem,1),
                b(Elem(j,:))=b(Elem(j,:))+det([1 1 1; Coord(Elem(j,:),:)'])*...
                    f([sum(Coord(Elem(j,:),:))/3,k*dt])/6;
            end   
            b=b(FreeNodes);
            
            yh(FreeNodes)=(1.5*M_h-dt*(A_h))\(M_h*(2*y_temp1-0.5*y_temp)+dt*b);
            y_temp=y_temp1;
            y_temp1=yh(FreeNodes);  
        end
    
    if (nl > 1)
        [L2e(nl-1),H1e(nl-1)]= ...
        Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempy, yh);
    end
        
        tempy=yh;
        Coord1=Coord;
        Elem1=Elem;
        
%         z=z_nodes(Coord,T); 
%         [L2e(nl),H1e(nl)]=Err(Coord, Elem, yh, z,T);
    
 
end

for ml=1:(nl-2)
    ocl2(ml+1)=abs(log(L2e(ml)/L2e(ml+1))/log(h(ml)/h(ml+1)));
    och1(ml+1)=abs(log(H1e(ml)/H1e(ml+1))/log(h(ml)/h(ml+1)));
end
ocl2(1)  = NaN;
och1(1)  = NaN;


format long e
Order_of_Convergence=table(L2e', ocl2', H1e', och1');
Order_of_Convergence.Properties.VariableNames = {'L2_Err_y' 'L2_Conv_y'...
                'H1_Err_y' 'H1_Conv_y' }
            
figure;
loglog(h(2:end), L2e, 'sb-', 'LineWidth', 1.1, 'MarkerSize', 8);
hold on;
loglog(h(2:end), H1e, 'sr-', 'LineWidth', 1.1, 'MarkerSize', 8);
grid on;

xlabel('Discretization parameter h');
ylabel('Errors');
title('Log–Log plot of errors aginst h');
legend('err_{L^2}(z_{h_i})', 'err_{H^1}(z_{h_i})');
