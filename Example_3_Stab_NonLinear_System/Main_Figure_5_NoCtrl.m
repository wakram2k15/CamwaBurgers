clear all;
format long;

% Geometry of finite elements - triangulation, local global node numbering,
% coordinates, boundary nodes

[Coord,Elem,Nb,Db]=InitialMesh(1);

T=1;
shift=25; eta=1;
for nl=1:3
   
    %% Edge-Node-Element Connections
    [n2ed,ed2el]=edge(Elem,Coord);
    %% Element Redrefine
    [Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);

    %%%%%%%%%%% find h %%%%%%%%%%%%%%
    h(nl)=sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
    dt=0.001;
 
    %% FullNodes and FreeNodes, notations of size
    FullNodes=[1:size(Coord,1)]; % Global node numbers
    FreeNodes=setdiff(FullNodes, unique(Db));
 
    nfull=length(FullNodes);
    nfree=length(FreeNodes); 
   
    M_1=sparse(size(Coord,1),size(Coord,1)); 
    K_1=sparse(size(Coord,1),size(Coord,1));
    A_1=sparse(size(Coord,1),size(Coord,1));
    A_2=sparse(size(Coord,1),size(Coord,1));

    for j=1:size(Elem,1)
        [ST,Mass,Ah1]=stima_1(Coord(Elem(j,:),:));
        K_1(Elem(j,:),Elem(j,:))=K_1(Elem(j,:),Elem(j,:))+ST;
        M_1(Elem(j,:),Elem(j,:))=M_1(Elem(j,:),Elem(j,:))+Mass; 
        A_1(Elem(j,:),Elem(j,:))=A_1(Elem(j,:),Elem(j,:))+ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
        [ysxev ysyev]=ysxe([sum(Coord(Elem(j,:),:))/3,0]);
        A_2(Elem(j,:),Elem(j,:))=A_2(Elem(j,:),Elem(j,:))+(ysxev+ysyev)*Mass;
    end
    
    Awh=-eta*K_1(FreeNodes,FreeNodes)-A_1(FreeNodes,FreeNodes)-A_2(FreeNodes,FreeNodes)+...
                shift*M_1(FreeNodes,FreeNodes);
    E1=M_1(FreeNodes,FreeNodes);


       
    %%%%%%%%%%%%%%%%%%%
    z_temp=zeros(length(FullNodes),1);
        for j=1:length(FullNodes)
           z_temp(FullNodes(j),1) = z_0([Coord(FullNodes(j),:),0]);
        end
   Ez(1,nl) =Norm(Coord, Elem, z_temp);
        
%   z_iguess(FreeNodes)=z_temp(FreeNodes); 


    for k=2  %% Time loop
         
        z_rand=rand(length(FullNodes),1);
        z_temp1=sparse(length(FullNodes),1);
        z_temp1(FreeNodes)=z_temp(FreeNodes);
        
        error= Norm(Coord, Elem, z_temp);
        while error>0.0001 %%Newton's loop
            
            M=sparse(size(Coord,1),size(Coord,1)); 
            K=sparse(size(Coord,1),size(Coord,1));
            A1=sparse(size(Coord,1),size(Coord,1));
            A2=sparse(size(Coord,1),size(Coord,1));
            Nz=sparse(size(Coord,1),size(Coord,1));
            Vz=sparse(size(Coord,1),size(Coord,1));

            for j=1:size(Elem,1)
                zt=z_temp1(Elem(j,:));
                [ST,Mass,Ah1,Nh,V]=stima_all(Coord(Elem(j,:),:),zt);
                K(Elem(j,:),Elem(j,:))=K(Elem(j,:),Elem(j,:))+ST;
                M(Elem(j,:),Elem(j,:))=M(Elem(j,:),Elem(j,:))+Mass;  
                A1(Elem(j,:),Elem(j,:))=A1(Elem(j,:),Elem(j,:))+ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                [ysxev ysyev]=ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                A2(Elem(j,:),Elem(j,:))=A2(Elem(j,:),Elem(j,:))+(ysxev+ysyev)*Mass;
                Nz(Elem(j,:),Elem(j,:))=Nz(Elem(j,:),Elem(j,:))+Nh;
                Vz(Elem(j,:),Elem(j,:))=Vz(Elem(j,:),Elem(j,:))+Nh+V;
            end
%     
            %% Construction of the matricse A, B where EY'=AY+Bu, Y=(y,z)^T at FullNodes
            Awh=-eta*K(FreeNodes,FreeNodes)-A1(FreeNodes,FreeNodes)-A2(FreeNodes,FreeNodes)+...
                shift*M(FreeNodes,FreeNodes);
            G =M(FreeNodes,FreeNodes)-dt*Awh-dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes);
            F = -( M(FreeNodes,FreeNodes)-dt*Awh...
                - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes))*z_temp1(FreeNodes)...
                + M(FreeNodes,FreeNodes)*z_temp(FreeNodes);
            zh=zeros(length(FullNodes),1);
            zh(FreeNodes)=G\F;
    
            error=Norm(Coord, Elem, zh);
            
            z_temp1(FreeNodes)=z_temp1(FreeNodes)+zh(FreeNodes);
            
        end
        
        Ez(k,nl)=Norm(Coord, Elem, z_temp1);
        
    end


    for k=3:T/dt+1  %% Time loop BDF2
         
        z_rand=rand(length(FullNodes),1);
        z_iguess=sparse(length(FullNodes),1);
        z_iguess(FreeNodes)=z_temp1(FreeNodes);
        
        error= Norm(Coord, Elem, z_temp);
        while error>0.0001 %%Newton's loop
            
            M=sparse(size(Coord,1),size(Coord,1)); 
            K=sparse(size(Coord,1),size(Coord,1));
            A1=sparse(size(Coord,1),size(Coord,1));
            A2=sparse(size(Coord,1),size(Coord,1));
            Nz=sparse(size(Coord,1),size(Coord,1));
            Vz=sparse(size(Coord,1),size(Coord,1));

            for j=1:size(Elem,1)
                zt=z_iguess(Elem(j,:));
                [ST,Mass,Ah1,Nh,V]=stima_all(Coord(Elem(j,:),:),zt);
                K(Elem(j,:),Elem(j,:))=K(Elem(j,:),Elem(j,:))+ST;
                M(Elem(j,:),Elem(j,:))=M(Elem(j,:),Elem(j,:))+Mass;  
                A1(Elem(j,:),Elem(j,:))=A1(Elem(j,:),Elem(j,:))+ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                [ysxev ysyev]=ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                A2(Elem(j,:),Elem(j,:))=A2(Elem(j,:),Elem(j,:))+(ysxev+ysyev)*Mass;
                Nz(Elem(j,:),Elem(j,:))=Nz(Elem(j,:),Elem(j,:))+Nh;
                Vz(Elem(j,:),Elem(j,:))=Vz(Elem(j,:),Elem(j,:))+Nh+V;
            end
%     
            %% Construction of the matricse A, B where EY'=AY+Bu, Y=(y,z)^T at FullNodes
            Awh=-eta*K(FreeNodes,FreeNodes)-A1(FreeNodes,FreeNodes)-A2(FreeNodes,FreeNodes)+...
                shift*M(FreeNodes,FreeNodes);
            G = M(FreeNodes,FreeNodes)-dt*Awh-dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes);
            F = -( 1.5*M(FreeNodes,FreeNodes)-dt*Awh...
                - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes))*z_iguess(FreeNodes)...
                +2*M(FreeNodes,FreeNodes)*z_temp1(FreeNodes)-0.5*M(FreeNodes,FreeNodes)*z_temp(FreeNodes);
            zh=zeros(length(FullNodes),1);
            zh(FreeNodes)=G\F;
    
            error=Norm(Coord, Elem, zh);
            
            z_iguess(FreeNodes)=z_iguess(FreeNodes)+zh(FreeNodes);
            
        end
        z_temp(FreeNodes)=z_temp1(FreeNodes);
        z_temp1(FreeNodes)=z_iguess(FreeNodes);
        
        Ez(k,nl)=Norm(Coord, Elem, z_iguess);
        
    end
  %  z_temp(FreeNodes)=z_iguess(FreeNodes);
    
        if (nl > 1)
            [L2ez(nl-1),H1ez(nl-1)]= Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempz, z_iguess);
        end
        
        tempz=z_iguess;
        Coord1=Coord;
        Elem1=Elem; 
        
 end
% 

%%% Energy plot
figure(1)
Time=0.01:0.01:T;
ee=0.01/dt:0.01/dt:T/dt;
Time=0:0.001:T;
ee=1:1:T/dt+1;
plot(Time, Ez(ee,1), 'g-'); hold on
plot(Time, Ez(ee,2), 'r-'); 
plot(Time, Ez(ee,3), 'b-');
% plot(Time, Ez(ee,4), 'c-');
% plot(Time, Ez(ee,5), 'k-');
xlabel('Time $t$'); ylabel('$\|z_h(t)\|$', 'Interpreter', 'latex');
legend('h=2^{-2}','h=2^{-3}','h=2^{-4}');
title('Evolution of solution (without control) with time $t$', 'Interpreter', 'latex');
grid on;

%% Semilog plotting - comment Energy Plotting
figure(2)
Time=0.01:0.01:T;
ee=0.01/dt:0.01/dt:T/dt;
Time=0:0.001:T;
ee=1:1:T/dt+1;
semilogy(Time, Ez(:,1)/Ez(1,1), '-g'); hold on
semilogy(Time, Ez(:,2)/Ez(1,2), '-r');
semilogy(Time, Ez(:,3)/Ez(1,3), '-b');
% semilogy(Time, Ez(:,4)/Ez(1,4), '-c');
% semilogy(Time, Ez(:,5)/Ez(1,5), '-m');
grid on;
title('Evolution of solution (without control) on log-scale with time t', 'Interpreter', 'latex')
xlabel('Time $t$', 'Interpreter', 'latex')
ylabel('Energy on log-scale', 'Interpreter', 'latex')
