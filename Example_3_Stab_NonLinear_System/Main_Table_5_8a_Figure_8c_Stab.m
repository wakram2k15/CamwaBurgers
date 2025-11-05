clear all;
format long;

[Coord,Elem,Nb,Db]=InitialMesh(1);

T=0.1;
shift=25; eta=1;
total_start = tic;
for nl=1:4
   refinement_start = tic; 
    %% Edge-Node-Element Connections
    [n2ed,ed2el]=edge(Elem,Coord);
    %% Element Redrefine
    [Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);

    %%%%%%%%%%% find h %%%%%%%%%%%%%%
    h(nl)=sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
    dt=0.001;
 
%   %% FullNodes and FreeNodes, notations of size
    FullNodes=[1:size(Coord,1)]; % Global node numbers
    FreeNodes=setdiff(FullNodes, unique(Db));
 
    nfull=length(FullNodes);
    nfree=length(FreeNodes); 
   
%%%%%Feedback Matricx calculation
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


%% For Stabilization

B=E1;
ne = nl;      % how many eigenvalues to compute
nc = 1;

%opts.p = 4*ne; % number of Lanczos vectors
 [V,D] = eigs(Awh,E1, nfree-2, 'SM'); % compute eigenvalues of (A,M)
 [Z,Dt] = eigs(Awh',E1, nfree-2, 'SM'); % compute eigenvalues of (A,M)
%save('eigShift6.mat','V','D')
%disp('Eigenvalues of A1')
D = diag(D);
Dt = diag(Dt);
% Z = V;
% Zt = Vt;
[d,ii]=sort(real(D), 'descend');
[dp,iit]=sort(real(Dt), 'descend');
%disp('20 eigenvalues of A1 and A1^T with largest real part')
D1 = D(ii(1:ne));
Dt = D(ii(1:ne));
% find unstable eig
iu = find(real(D1) > 0);
nu = length(iu);
fprintf(1, 'Number of unstable eigenvalues of A = %d\n', nu)

% with shift > 5 first pair of D1 must be unstable
assert(min(real(D(ii(1:nu)))) > 0.0)
Vt = V(:,ii(1:ne));
Zt = Z(:,iit(1:ne));

% Make Vt and Zt orthonormal
% p must be diagonal
p = Vt.' * E1 * Zt;
% Check p is diagonal and diagonal entries are non-zero
assert(min(abs(diag(p))) > 0.0);
%assert(is_diag(p)==1)
p = diag(p);

% normalize
for j=1:ne
   Zt(:,j) = Zt(:,j) / p(j);
end
% 
Vy = Vt(:,1:nu);
Zy = Zt(:,1:nu);


Au = Zy' * Awh * Vy;
Bu = Zy'*B;
Qu = Vy'*E1*Vy;
Ru = eye(nfree);
ric_start = tic;
[Pu,L,G]=care(Au,Bu,Qu,Ru); % Solves Riccati equation
time_riccati(nl) = toc(ric_start);

% Feedback operator
FB = Ru \ ((B' * Zy) * Pu * (Zy' * E1));
%FB=B'*Pu;
[Vs,Ds] = eigs(Awh-B*FB,E1, nfree-2, 'SM');
Ds = diag(Ds);
[ds,iis]=sort(real(Ds), 'descend');
       
    %%%%%%%%%%%%%%%%%%%
    z_temp=zeros(length(FullNodes),1);
        for j=1:length(FullNodes)
           z_temp(FullNodes(j),1) = z_0([Coord(FullNodes(j),:),0]);
        end
   Ez(1,nl) =Norm(Coord, Elem, z_temp);
        
%   z_iguess(FreeNodes)=z_temp(FreeNodes); 
   newton_time = 0;
    for k=2:T/dt  %% Time loop
         
        z_rand=rand(length(FullNodes),1);
        z_iguess=sparse(length(FullNodes),1);
        z_iguess(FreeNodes)=z_temp(FreeNodes);
        
        error= Norm(Coord, Elem, z_temp);
        newton_iter_start = tic;
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
            G = M(FreeNodes,FreeNodes)-dt*Awh-dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes)+dt*B*FB;
            F = -( M(FreeNodes,FreeNodes)-dt*Awh...
                - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes)+dt*B*FB)*z_iguess(FreeNodes)...
                + M(FreeNodes,FreeNodes)*z_temp(FreeNodes);
            zh=zeros(length(FullNodes),1);
            zh(FreeNodes)=G\F;
    
            error=Norm(Coord, Elem, zh);
            
            z_iguess(FreeNodes)=z_iguess(FreeNodes)+zh(FreeNodes);
            
        end
        z_temp(FreeNodes)=z_iguess(FreeNodes);
        newton_time = newton_time + toc(newton_iter_start);
        Ez(k,nl)=Norm(Coord, Elem, z_iguess);
        uh=zeros(length(FullNodes),1);
        uh(FreeNodes)=-FB*z_iguess(FreeNodes);
        Eu(k,nl)=Norm(Coord, Elem, uh);
        
    end
  %  z_temp(FreeNodes)=z_iguess(FreeNodes);
  time_newton(nl) = newton_time;
    
        if (nl > 1)
            [L2ez(nl-1),H1ez(nl-1)]= Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempz, z_iguess);
            [L2eu(nl-1),H1eu(nl-1)]= Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempu, uh);
        end
        
        tempz=z_iguess;
        tempu=uh;
        Coord1=Coord;
        Elem1=Elem; 
        
        time_refine(nl) = toc(refinement_start);
 end
% 
total_time = toc(total_start);
%  
for pl=1:(nl-2)
     oczl2(pl+1)=log(L2ez(pl)/L2ez(pl+1))/log(h(pl)/h(pl+1));
     oczh1(pl+1)=log(H1ez(pl)/H1ez(pl+1))/log(h(pl)/h(pl+1));
     ocul2(pl+1)=log(L2eu(pl)/L2eu(pl+1))/log(h(pl)/h(pl+1));
end

oczl2(1)  = NaN;
oczh1(1)  = NaN;
ocul2(1) = NaN;

format long e
Order_of_Convergence=table(L2ez', oczl2', H1ez', oczh1', L2eu', ocul2');
Order_of_Convergence.Properties.VariableNames = {'L2_Err_y' 'L2_Conv_y'...
                'H1_Err_y' 'H1_Conv_y' 'L2_Err_u' 'L2_Conv_u' }


%%%Log-log plot
loglog(h(2:end), L2ez, 'sb-', 'LineWidth', 1.1, 'MarkerSize', 8);
hold on;
loglog(h(2:end), L2eu, 'sk-', 'LineWidth', 1.1, 'MarkerSize', 8);
loglog(h(2:end), H1ez, 'sr-', 'LineWidth', 1.1, 'MarkerSize', 8);
grid on;

xlabel('Discretization parameter h');
ylabel('Errors');
title('Log–Log plot of errors aginst discretization parameter h');
legend('err_{L^2}(z_{h_i})', 'err_{L^2}(u_{h_i})', 'err_{H^1}(z_{h_i})');


disp('Timing summary:');
for nl=1:4
    fprintf('Refinement %d:\n', nl);
    %fprintf('   Eigenvalue computation: %.4f sec\n', time_eigs(nl));
    fprintf('   Riccati solver:         %.4f sec\n', time_riccati(nl));
    fprintf('   Newton iterations:      %.4f sec\n', time_newton(nl));
    fprintf('   TOTAL per refinement:   %.4f sec\n\n', time_refine(nl));
end
fprintf('TOTAL run time: %.4f sec\n', total_time);
