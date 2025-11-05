% Projection-based (modal) stabilization for viscous Burgers (adapted to your code)
clear all; close all; format long;

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
    dt=0.01;   % chosen as in your working runs

    FullNodes=[1:size(Coord,1)]; % Global node numbers
    FreeNodes=setdiff(FullNodes, unique(Db));

    nfull = length(FullNodes);
    nfree = length(FreeNodes);
    
   %%%%Construction of FE matrices  
    M_1 = sparse(size(Coord,1),size(Coord,1)); 
    K_1 = sparse(size(Coord,1),size(Coord,1));
    A_1 = sparse(size(Coord,1),size(Coord,1));
    A_2 = sparse(size(Coord,1),size(Coord,1));

    for j=1:size(Elem,1)
        [ST,Mass,Ah1]=stima_1(Coord(Elem(j,:),:));
        K_1(Elem(j,:),Elem(j,:)) = K_1(Elem(j,:),Elem(j,:)) + ST;
        M_1(Elem(j,:),Elem(j,:)) = M_1(Elem(j,:),Elem(j,:)) + Mass; 
        A_1(Elem(j,:),Elem(j,:)) = A_1(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
        [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
        A_2(Elem(j,:),Elem(j,:)) = A_2(Elem(j,:),Elem(j,:)) + (ysxev+ysyev)*Mass;
    end
    
    % linearized operator on free DOFs
    Awh = -eta*K_1(FreeNodes,FreeNodes) - A_1(FreeNodes,FreeNodes) - A_2(FreeNodes,FreeNodes) + ...
                shift*M_1(FreeNodes,FreeNodes);
    E1 = M_1(FreeNodes,FreeNodes);

    % ========== actuator / control operator B ==========
    B = E1;      

    % ===================== PROJECTION / MODAL REDUCTION =====================
    ne_compute = min(80, nfree);    % cap to avoid asking too many eigenpairs
    opts.tol = 1e-6;
    try
        [Vfull, Dfull] = eigs(Awh, E1, ne_compute-2, 'SM', opts);  % nfree x ne_compute
        [Zfull, Dfull2] = eigs(Awh', E1, ne_compute-2, 'SM', opts); % adjoint
    catch
        % fallback to dense eig if eigs fails (small nfree)
        [Vfull_all, D_all] = eig(Awh, E1);
        [Zfull_all, D2_all] = eig(Awh', E1);
        Vfull = Vfull_all;
        Zfull = Zfull_all;
    end
    

    % sort by real part descending
    dvals = diag( (Vfull' * (Awh * Vfull)) ./ (Vfull' * (E1 * Vfull) + eps ) ); 
    [~, idx_sort] = sort(real(dvals), 'descend');
    Vt = Vfull(:, idx_sort);
    Zt = Zfull(:, idx_sort);

    ne_use = min( min(20, size(Vt,2)), size(Vt,2) );  
    Vt = Vt(:,1:ne_use);
    Zt = Zt(:,1:ne_use);

    pmat = Vt' * (E1 * Zt);
    pdiag = diag(pmat);
    smallp = find(abs(pdiag) < 1e-12);
    if ~isempty(smallp)
        warning('Small diagonal entries in pairing p; modes may be nearly orthogonal.');
        pdiag(smallp) = pdiag(smallp) + 1e-12;
    end
    for j=1:ne_use
        Zt(:,j) = Zt(:,j) / pdiag(j);
    end

    lambdas = zeros(ne_use,1);
    for j=1:ne_use
        lambdas(j) = (Zt(:,j)' * (Awh * Vt(:,j))) / (Zt(:,j)' * (E1 * Vt(:,j)));
    end
    unstable_idx = find(real(lambdas) > 1e-12);
    nu = length(unstable_idx);
    if nu == 0
        nu = min(2, ne_use);
        unstable_idx = 1:nu;
    end
    Vy = Vt(:, unstable_idx);   
    Zy = Zt(:, unstable_idx);   

    Au = Zy' * (Awh * Vy);      
    Bu = Zy' * B;               
    Qu = Vy' * (E1 * Vy);       

    m_act = size(B,2);
    Ru_reduced = eye(m_act);    
    % ===================== Riccati Solver =====================
    ric_start = tic;
    try
        [Pu,~,~] = care(Au, Bu, Qu, Ru_reduced);
    catch ME
        warning('Reduced CARE failed: %s. Trying regularization.', ME.message);
        Pu = care(Au, Bu, Qu + 1e-8*eye(size(Qu)), Ru_reduced);
    end

    FB = Ru_reduced \ ( (B' * Zy) * Pu * (Zy' * E1) );   
    BFB = B * FB;   

    % ==================== INITIAL DATA PROJECTION ====================
    z_temp=zeros(length(FullNodes),1);
    for j=1:size(Elem,1)
        z_temp(Elem(j,:)) = z_temp(Elem(j,:)) + det([1 1 1; Coord(Elem(j,:),:)'])*...
                    z_0([sum(Coord(Elem(j,:),:))/3,0])/6;
    end
    z_temp(FreeNodes) = E1 \ z_temp(FreeNodes);
    Ez(1,nl) = Norm(Coord, Elem, z_temp);

    % ==================== FIRST STEP (k = 2) ====================
    
    for k=2
        z_temp1 = sparse(length(FullNodes),1);
        z_temp1(FreeNodes) = z_temp(FreeNodes);

        error = Norm(Coord, Elem, z_temp);
        while error>1.0e-4
            M = sparse(size(Coord,1),size(Coord,1)); 
            K = sparse(size(Coord,1),size(Coord,1));
            A1loc = sparse(size(Coord,1),size(Coord,1));
            A2loc = sparse(size(Coord,1),size(Coord,1));
            Nz = sparse(size(Coord,1),size(Coord,1));
            Vz = sparse(size(Coord,1),size(Coord,1));

            for j=1:size(Elem,1)
                zt = z_temp1(Elem(j,:));
                [ST,Mass,Ah1,Nh,Vv] = stima_all(Coord(Elem(j,:),:),zt);
                K(Elem(j,:),Elem(j,:)) = K(Elem(j,:),Elem(j,:)) + ST;
                M(Elem(j,:),Elem(j,:)) = M(Elem(j,:),Elem(j,:)) + Mass;  
                A1loc(Elem(j,:),Elem(j,:)) = A1loc(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                A2loc(Elem(j,:),Elem(j,:)) = A2loc(Elem(j,:),Elem(j,:)) + (ysxev+ysyev)*Mass;
                Nz(Elem(j,:),Elem(j,:)) = Nz(Elem(j,:),Elem(j,:)) + Nh;
                Vz(Elem(j,:),Elem(j,:)) = Vz(Elem(j,:),Elem(j,:)) + Nh + Vv;
            end

            Awh_loc = -eta*K(FreeNodes,FreeNodes) - A1loc(FreeNodes,FreeNodes) - A2loc(FreeNodes,FreeNodes) + ...
                shift*M(FreeNodes,FreeNodes);

            G = M(FreeNodes,FreeNodes) - dt*Awh_loc - dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes) + dt*BFB;
            F = -( M(FreeNodes,FreeNodes) - dt*Awh_loc - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes) + dt*BFB ) * z_temp1(FreeNodes) ...
                + M(FreeNodes,FreeNodes)*z_temp(FreeNodes);

            zh = zeros(length(FullNodes),1);
            zh(FreeNodes) = G \ F;

            error = Norm(Coord, Elem, zh);
            z_temp1(FreeNodes) = z_temp1(FreeNodes) + zh(FreeNodes);
        end

        Ez(k,nl) = Norm(Coord, Elem, z_temp1);

        % compute control coefficients and lift to nodal control field
        u_coeff = -FB * z_temp1(FreeNodes);   
        uh = zeros(length(FullNodes),1);
        uh(FreeNodes) = B * u_coeff;         
        Eu(k,nl) = Norm(Coord, Elem, uh);
    end

    % ==================== BDF2 STEPS (k >= 3) ====================
    for k=3:T/dt
        z_iguess = sparse(length(FullNodes),1);
        z_iguess(FreeNodes) = z_temp1(FreeNodes);

        error = Norm(Coord, Elem, z_temp);
        newton_iter_start = tic;
        while error>1.0e-4
            M = sparse(size(Coord,1),size(Coord,1)); 
            K = sparse(size(Coord,1),size(Coord,1));
            A1loc = sparse(size(Coord,1),size(Coord,1));
            A2loc = sparse(size(Coord,1),size(Coord,1));
            Nz = sparse(size(Coord,1),size(Coord,1));
            Vz = sparse(size(Coord,1),size(Coord,1));

            for j=1:size(Elem,1)
                zt = z_iguess(Elem(j,:));
                [ST,Mass,Ah1,Nh,Vv] = stima_all(Coord(Elem(j,:),:),zt);
                K(Elem(j,:),Elem(j,:))  = K(Elem(j,:),Elem(j,:))  + ST;
                M(Elem(j,:),Elem(j,:))  = M(Elem(j,:),Elem(j,:))  + Mass;  
                A1loc(Elem(j,:),Elem(j,:)) = A1loc(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                A2loc(Elem(j,:),Elem(j,:)) = A2loc(Elem(j,:),Elem(j,:)) + (ysxev+ysyev)*Mass;
                Nz(Elem(j,:),Elem(j,:)) = Nz(Elem(j,:),Elem(j,:)) + Nh;
                Vz(Elem(j,:),Elem(j,:)) = Vz(Elem(j,:),Elem(j,:)) + Nh + Vv;
            end

            Awh_loc = -eta*K(FreeNodes,FreeNodes) - A1loc(FreeNodes,FreeNodes) - A2loc(FreeNodes,FreeNodes) + ...
                shift*M(FreeNodes,FreeNodes);

            G = M(FreeNodes,FreeNodes) - dt*Awh_loc - dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes) + dt*BFB;
            F = -( 1.5*M(FreeNodes,FreeNodes) - dt*Awh_loc - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes) + dt*BFB ) * z_iguess(FreeNodes) ...
                + 2*M(FreeNodes,FreeNodes)*z_temp1(FreeNodes) - 0.5*M(FreeNodes,FreeNodes)*z_temp(FreeNodes);

            zh = zeros(length(FullNodes),1);
            zh(FreeNodes) = G \ F;

            error = Norm(Coord, Elem, zh);
            z_iguess(FreeNodes) = z_iguess(FreeNodes) + zh(FreeNodes);
        end
        

        z_temp(FreeNodes) = z_temp1(FreeNodes);
        z_temp1(FreeNodes)= z_iguess(FreeNodes);

        Ez(k,nl) = Norm(Coord, Elem, z_iguess);

        % compute control coefficients and lift
        u_coeff = -FB * z_temp(FreeNodes);
        uh = zeros(length(FullNodes),1);
        uh(FreeNodes) = B * u_coeff;
        Eu(k,nl) = Norm(Coord, Elem, uh);
    end

    % convergence checks (unchanged)
    if (nl > 1)
        [L2ez(nl-1),H1ez(nl-1)] = Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempz, z_iguess);
        [L2eu(nl-1),H1eu(nl-1)] = Err_gaussiantest(Coord,Coord1, Elem1,ed2el,n2ed ,tempu, uh);
    end
    tempz = z_iguess;
    tempu = uh;
    Coord1 = Coord;
    Elem1 = Elem;

end


%%% Energy Plotting
      figure(1);
        Time=dt:dt:T;
        plot(Time, Ez(:,1), '-g'); hold on
        plot(Time, Ez(:,2), '-r');
        plot(Time, Ez(:,3), '-b');
%        plot(Time, Ez(:,4), '-c');
%         plot(Time, Ez(:,5), '-m');
%         plot(Time, Ez(:,6), '-k');
        grid on;
        title('Evolution of stabilized solution against time $t$',...
            'Interpreter', 'latex')
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\|z_h(t)\|$', 'Interpreter', 'latex')
        
%%% Control Norm Plotting
        figure(2)
        Time=dt:dt:T;
        plot(Time, Eu(:,1), '-g'); hold on
        plot(Time, Eu(:,2), '-r');
        plot(Time, Eu(:,3), '-b');
%        plot(Time, Eu(:,4), '-c');
%         plot(Time, Eu(:,5), '-m');
%         plot(Time, Eu(:,6), '-k');
        grid on;
        title('Evolution of stabilizing control against time $t$', 'Interpreter', 'latex')
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\|u_h(t)\|$', 'Interpreter', 'latex')




