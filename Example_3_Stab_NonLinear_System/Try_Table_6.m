clear all;
format long;

% Geometry of finite elements - triangulation, local/global node numbering,
% coordinates, boundary nodes
[Coord,Elem,Nb,Db] = InitialMesh(1);

T = 1;
shift = 25; 
eta   = 1;

for nl = 1:3
   
    %% Edge-Node-Element Connections
    [n2ed,ed2el] = edge(Elem,Coord);

    %% Element Red-refine
    [Coord,Elem,Db,Nb] = redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);

    % -------- mesh size & time step --------
    h(nl) = sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
    dt    = 0.01;  % or dt = 0.49*h(nl)^2;

    % -------- nodes --------
    FullNodes = 1:size(Coord,1);            % global node numbers
    FreeNodes = setdiff(FullNodes, unique(Db));
    nfull     = length(FullNodes);
    nfree     = length(FreeNodes);
    
    % -------- FE matrices for feedback construction --------
    M_1 = sparse(size(Coord,1),size(Coord,1)); 
    K_1 = sparse(size(Coord,1),size(Coord,1));
    A_1 = sparse(size(Coord,1),size(Coord,1));
    A_2 = sparse(size(Coord,1),size(Coord,1));

    for j = 1:size(Elem,1)
        [ST,Mass,Ah1] = stima_1(Coord(Elem(j,:),:));
        K_1(Elem(j,:),Elem(j,:)) = K_1(Elem(j,:),Elem(j,:)) + ST;
        M_1(Elem(j,:),Elem(j,:)) = M_1(Elem(j,:),Elem(j,:)) + Mass; 
        A_1(Elem(j,:),Elem(j,:)) = A_1(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
        [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
        A_2(Elem(j,:),Elem(j,:)) = A_2(Elem(j,:),Elem(j,:)) + (ysxev+ysyev)*Mass;
    end
    
    Awh = -eta*K_1(FreeNodes,FreeNodes) ...
          -A_1(FreeNodes,FreeNodes) ...
          -A_2(FreeNodes,FreeNodes) ...
          + shift*M_1(FreeNodes,FreeNodes);

    E1  = M_1(FreeNodes,FreeNodes);

    % ==================== LOCALIZED INTERIOR CONTROL ====================
    % Assemble mass matrix only over elements fully inside omega
    M_loc = sparse(size(Coord,1),size(Coord,1));
    for j = 1:size(Elem,1)
        centroid = mean(Coord(Elem(j,:),:));
        if centroid(1) >= 0.10 && centroid(1) <= 0.90 && ...
           centroid(2) >= 0.10 && centroid(2) <= 0.90
            [~,Mass] = stima_1(Coord(Elem(j,:),:));
            M_loc(Elem(j,:),Elem(j,:)) = M_loc(Elem(j,:),Elem(j,:)) + Mass;
        end
    end
    % Restrict to free nodes
    B = M_loc(FreeNodes,FreeNodes);

    % CARE (use R = I; capture only Pu to avoid name clashes)
    Ru = eye(size(B,2));
    [Pu, ~, ~] = care((E1\Awh), (E1\B), E1, Ru);

    % Feedback operator FB
    FB = Ru \ ((E1\B)' * Pu);

    % ==================== INITIAL DATA PROJECTION ====================
    z_temp = zeros(length(FullNodes),1);
    for j = 1:size(Elem,1)
        z_temp(Elem(j,:)) = z_temp(Elem(j,:)) + ...
            det([1 1 1; Coord(Elem(j,:),:)']) * z_0([sum(Coord(Elem(j,:),:))/3,0]) / 6;
    end 
    z_temp(FreeNodes) = E1 \ z_temp(FreeNodes);
    Ez(1,nl)          = Norm(Coord, Elem, z_temp);
        
    % ==================== FIRST STEP (k = 2) ====================
    for k = 2  % Time loop (first step)
        z_temp1 = sparse(length(FullNodes),1);
        z_temp1(FreeNodes) = z_temp(FreeNodes);
        
        error = Norm(Coord, Elem, z_temp);
        while error > 1.0e-4  % Newton loop
            M  = sparse(size(Coord,1),size(Coord,1)); 
            K  = sparse(size(Coord,1),size(Coord,1));
            A1 = sparse(size(Coord,1),size(Coord,1));
            A2 = sparse(size(Coord,1),size(Coord,1));
            Nz = sparse(size(Coord,1),size(Coord,1));
            Vz = sparse(size(Coord,1),size(Coord,1));

            for j = 1:size(Elem,1)
                zt = z_temp1(Elem(j,:));
                [ST,Mass,Ah1,Nh,V] = stima_all(Coord(Elem(j,:),:),zt);
                K(Elem(j,:),Elem(j,:))  = K(Elem(j,:),Elem(j,:))  + ST;
                M(Elem(j,:),Elem(j,:))  = M(Elem(j,:),Elem(j,:))  + Mass;  
                A1(Elem(j,:),Elem(j,:)) = A1(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                A2(Elem(j,:),Elem(j,:)) = A2(Elem(j,:),Elem(j,:)) + (ysxev+ysyev)*Mass;
                Nz(Elem(j,:),Elem(j,:)) = Nz(Elem(j,:),Elem(j,:)) + Nh;
                Vz(Elem(j,:),Elem(j,:)) = Vz(Elem(j,:),Elem(j,:)) + Nh + V;
            end

            Awh = -eta*K(FreeNodes,FreeNodes) ...
                  -A1(FreeNodes,FreeNodes) ...
                  -A2(FreeNodes,FreeNodes) ...
                  + shift*M(FreeNodes,FreeNodes);

            Gmat =  M(FreeNodes,FreeNodes) ...
                    - dt*Awh ...
                    - dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes) ...
                    + dt*B*FB;

            Fvec = -( M(FreeNodes,FreeNodes) ...
                      - dt*Awh ...
                      - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes) ...
                      + dt*B*FB ) * z_temp1(FreeNodes) ...
                   + M(FreeNodes,FreeNodes)*z_temp(FreeNodes);

            zh = zeros(length(FullNodes),1);
            zh(FreeNodes) = Gmat \ Fvec;
    
            error = Norm(Coord, Elem, zh);
            z_temp1(FreeNodes) = z_temp1(FreeNodes) + zh(FreeNodes);
        end
        
        Ez(k,nl) = Norm(Coord, Elem, z_temp1);

        % Localized control (variational)
        u_ctrl = -FB * z_temp1(FreeNodes);
        uh     = zeros(length(FullNodes),1);
        uh(FreeNodes) = u_ctrl;  % variational lifting
        Eu(k,nl) = Norm(Coord, Elem, uh);
    end

    % ==================== BDF2 STEPS (k >= 3) ====================
    for k = 3:T/dt+1
        z_iguess = sparse(length(FullNodes),1);
        z_iguess(FreeNodes) = z_temp1(FreeNodes);
        
        error = Norm(Coord, Elem, z_temp);
        while error > 1.0e-4  % Newton loop
            M  = sparse(size(Coord,1),size(Coord,1)); 
            K  = sparse(size(Coord,1),size(Coord,1));
            A1 = sparse(size(Coord,1),size(Coord,1));
            A2 = sparse(size(Coord,1),size(Coord,1));
            Nz = sparse(size(Coord,1),size(Coord,1));
            Vz = sparse(size(Coord,1),size(Coord,1));

            for j = 1:size(Elem,1)
                zt = z_iguess(Elem(j,:));
                [ST,Mass,Ah1,Nh,V] = stima_all(Coord(Elem(j,:),:),zt);
                K(Elem(j,:),Elem(j,:))  = K(Elem(j,:),Elem(j,:))  + ST;
                M(Elem(j,:),Elem(j,:))  = M(Elem(j,:),Elem(j,:))  + Mass;  
                A1(Elem(j,:),Elem(j,:)) = A1(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                A2(Elem(j,:),Elem(j,:)) = A2(Elem(j,:),Elem(j,:)) + (ysxev+ysyev)*Mass;
                Nz(Elem(j,:),Elem(j,:)) = Nz(Elem(j,:),Elem(j,:)) + Nh;
                Vz(Elem(j,:),Elem(j,:)) = Vz(Elem(j,:),Elem(j,:)) + Nh + V;
            end

            Awh = -eta*K(FreeNodes,FreeNodes) ...
                  -A1(FreeNodes,FreeNodes) ...
                  -A2(FreeNodes,FreeNodes) ...
                  + shift*M(FreeNodes,FreeNodes);

            Gmat =  M(FreeNodes,FreeNodes) ...
                    - dt*Awh - dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes) ...
                    + dt*B*FB;

            Fvec = -( 1.5*M(FreeNodes,FreeNodes) ...
                      - dt*Awh - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes) ...
                      + dt*B*FB ) * z_iguess(FreeNodes) ...
                   + 2*M(FreeNodes,FreeNodes)*z_temp1(FreeNodes) ...
                   - 0.5*M(FreeNodes,FreeNodes)*z_temp(FreeNodes);

            zh = zeros(length(FullNodes),1);
            zh(FreeNodes) = Gmat \ Fvec;
    
            error = Norm(Coord, Elem, zh);
            z_iguess(FreeNodes) = z_iguess(FreeNodes) + zh(FreeNodes);
        end

        z_temp(FreeNodes)  = z_temp1(FreeNodes);
        z_temp1(FreeNodes) = z_iguess(FreeNodes);
        
        Ez(k,nl) = Norm(Coord, Elem, z_iguess);

        % Localized control (variational lifting)
        u_ctrl = -FB * z_temp(FreeNodes);
        uh     = zeros(length(FullNodes),1);
        uh(FreeNodes) = u_ctrl;
        Eu(k,nl) = Norm(Coord, Elem, uh);
    end
    
    % ===== convergence checks =====
    if (nl > 1)
        [L2ez(nl-1),H1ez(nl-1)] = Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempz, z_iguess);
        [L2eu(nl-1),H1eu(nl-1)] = Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempu, uh);
    end
    tempz = z_iguess;
    tempu = uh;
    Coord1 = Coord;
    Elem1  = Elem;
end

     %%% Energy Plotting
      figure(1);
        Time=dt:dt:T;
        plot(Time, Ez(:,1), '-g'); hold on
        plot(Time, Ez(:,2), '-r');
        plot(Time, Ez(:,3), '-b');
%         plot(Time, Ez(:,4), '-c');
%         plot(Time, Ez(:,5), '-m');
%         plot(Time, Ez(:,6), '-k');
        grid on;
        title('Evolution of stabilized solution $z_h^\sharp$ with time $t$',...
            'Interpreter', 'latex')
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\|z_h^\sharp(t)\|$', 'Interpreter', 'latex')
        
     %%% Semilog plotting - comment Energy Plotting
     figure(2);
        Time=dt:dt:T;
        semilogy(Time, Ez(:,1)/Ez(1,1), '-g'); hold on
        semilogy(Time, Ez(:,2)/Ez(1,2), '-r');
        semilogy(Time, Ez(:,3)/Ez(1,3), '-b');
%         semilogy(Time, Ez(:,4)/Ez(1,4), '-c');
%         semilogy(Time, Ez(:,5)/Ez(1,5), '-m');
%         semilogy(Time, Ez(:,6)/Ez(1,6), '-k');
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
%         plot(Time, Eu(:,4), '-c');
%         plot(Time, Eu(:,5), '-m');
%         plot(Time, Eu(:,6), '-k');
        grid on;
        title('Evolution of stabilizing control $u_h^\sharp$ with time $t$', 'Interpreter', 'latex')
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\|u_h(t)\|$', 'Interpreter', 'latex')

for pl = 1:(nl-2)
    oczl2(pl+1) = log(L2ez(pl)/L2ez(pl+1))/log(h(pl)/h(pl+1));
    oczh1(pl+1) = log(H1ez(pl)/H1ez(pl+1))/log(h(pl)/h(pl+1));
    ocul2(pl+1) = log(L2eu(pl)/L2eu(pl+1))/log(h(pl)/h(pl+1));
end

oczl2(1)  = NaN;
oczh1(1)  = NaN;
ocul2(1) = NaN;

format long e
Order_of_Convergence=table(L2ez', oczl2', H1ez', oczh1', L2eu', ocul2');
Order_of_Convergence.Properties.VariableNames = {'L2_Err_y' 'L2_Conv_y'...
                'H1_Err_y' 'H1_Conv_y' 'L2_Err_u' 'L2_Conv_u' }