clear all;
format long;

[Coord,Elem,Nb,Db]=InitialMesh(1);


T=1; shift=24; 
eta=1;

for nl=1:4
   
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
    
    All_inv=M_h\[A_h eye(nfree)];
    A11=M_h\A_h; B1=eye(nfree); 
    invRu=M_h\eye(nfree);
    Ru=eye(nfree)*M_h;
    [Pu,L,G]=care(A11,B1,M_h,Ru); %
    
    yh=zeros(length(FullNodes),1);
    y_temp=zeros(length(FullNodes),1);
        
    for j=1:size(Elem,1),
                y_temp(Elem(j,:))=y_temp(Elem(j,:))+...
                    det([1 1 1; Coord(Elem(j,:),:)'])*y_0([sum(Coord(Elem(j,:),:))/3,0])/6-...
                 det([1 1 1; Coord(Elem(j,:),:)'])*ys([sum(Coord(Elem(j,:),:))/3,0])/6;
    end 
        
    y_temp=M_h\y_temp(FreeNodes); 
    
%% For Stabilization

B=M_h;
 
% Feedback operator
FB=Ru\(B1'*Pu);

        uh=zeros(length(FullNodes),1);
        for k=1 %BDF1 in first time step
            y_temp1=(M_h-dt*A_h+dt*B*FB)\((M_h)*y_temp);
            yh(FreeNodes)=y_temp1;
            uh(FreeNodes)=-FB*y_temp1;
            E(k,nl)=Norm(Coord, Elem, yh);
            Eu(k,nl)=Norm(Coord, Elem, uh);
        end
%     
        for k=2:T/dt %BDF2 in second time step onwards
%                          
            yh(FreeNodes)=(1.5*M_h-dt*(A_h)+dt*B*FB)\(M_h*(2*y_temp1-0.5*y_temp));
            y_temp=y_temp1;
            y_temp1=yh(FreeNodes);
            %yh(FreeNodes)=yh(1:nfree,1);
            uh(FreeNodes)=-FB*y_temp1;
            E(k,nl)=Norm(Coord, Elem, yh);%+Norm(Coord, Elem, zh);    
            Eu(k,nl)=Norm(Coord, Elem, uh);
        end


    
   %  y=u_nodes(Coord,T);
    
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

      %%% Energy Plotting
      figure(1);
        Time=dt:dt:T;
        plot(Time, E(:,1), '-g'); hold on
        plot(Time, E(:,2), '-r');
        plot(Time, E(:,3), '-b');
        plot(Time, E(:,4), '-c');
%         plot(Time, E(:,5), '-m');
%         plot(Time, E(:,6), '-k');
        grid on;
        title('Evolution of stabilized solution $z_h^\sharp$ with time $t$',...
            'Interpreter', 'latex')
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\|z_h^\sharp(t)\|$', 'Interpreter', 'latex')
        
     %%% Semilog plotting - comment Energy Plotting
     figure(2);
        Time=dt:dt:T;
        semilogy(Time, E(:,1)/E(1,1), '-g'); hold on
        semilogy(Time, E(:,2)/E(1,2), '-r');
        semilogy(Time, E(:,3)/E(1,3), '-b');
        semilogy(Time, E(:,4)/E(1,4), '-c');
%         semilogy(Time, E(:,5)/E(1,5), '-m');
%         semilogy(Time, E(:,6)/E(1,6), '-k');
        grid on;
          title('Evolution of stabilized solution $z_h^\sharp$  on log-scale with time $t$',...
              'Interpreter', 'latex')
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\|z_h^\sharp(t)\|$ on log-scale', 'Interpreter', 'latex')


        %%% Control Norm Plotting
        figure(3)
        Time=dt:dt:T;
        plot(Time, Eu(:,1), '-g'); hold on
        plot(Time, Eu(:,2), '-r');
        plot(Time, Eu(:,3), '-b');
        plot(Time, Eu(:,4), '-c');
%         plot(Time, E(:,5), '-m');
%         plot(Time, E(:,6), '-k');
        grid on;
        title('Evolution of stabilizing control $u_h^\sharp$ with time $t$', 'Interpreter', 'latex')
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\|u_h(t)\|$', 'Interpreter', 'latex')







