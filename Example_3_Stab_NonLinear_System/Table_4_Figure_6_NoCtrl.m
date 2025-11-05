clear all;
format long;

%% =========================================================
%% MAIN DRIVER: runs both cases (T = 0.1 and T = 1)
%% =========================================================
for T = [0.1, 1]

    fprintf('\n==============================\n');
    fprintf(' Running simulation for T = %.1f\n', T);
    fprintf('==============================\n');

    shift = 25; eta = 1;

    % Geometry and initial mesh
    [Coord,Elem,Nb,Db] = InitialMesh(1);

    for nl = 1:3
        %% Edge-Node-Element Connections
        [n2ed,ed2el] = edge(Elem,Coord);
        %% Element Refine
        [Coord,Elem,Db,Nb] = redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);

        %%%%%%%%%%% find h %%%%%%%%%%%%%%
        h(nl) = sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
        dt = 0.001;
    
        %% FullNodes and FreeNodes
        FullNodes = 1:size(Coord,1);
        FreeNodes = setdiff(FullNodes, unique(Db));
    
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
        
        Awh = -eta*K_1(FreeNodes,FreeNodes) - A_1(FreeNodes,FreeNodes) - A_2(FreeNodes,FreeNodes) + ...
              shift*M_1(FreeNodes,FreeNodes);

        %% Initial condition
        z_temp = zeros(length(FullNodes),1);
        for j = 1:length(FullNodes)
           z_temp(FullNodes(j),1) = z_0([Coord(FullNodes(j),:),0]);
        end
        Ez(1,nl) = Norm(Coord, Elem, z_temp);
            
        for k = 2  %% Time loop
            z_temp1 = z_temp;
            error = Norm(Coord, Elem, z_temp);
            while error > 0.0001
                M = sparse(size(Coord,1),size(Coord,1)); 
                K = sparse(size(Coord,1),size(Coord,1));
                A1 = sparse(size(Coord,1),size(Coord,1));
                A2 = sparse(size(Coord,1),size(Coord,1));
                Nz = sparse(size(Coord,1),size(Coord,1));
                Vz = sparse(size(Coord,1),size(Coord,1));

                for j = 1:size(Elem,1)
                    zt = z_temp1(Elem(j,:));
                    [ST,Mass,Ah1,Nh,V] = stima_all(Coord(Elem(j,:),:),zt);
                    K(Elem(j,:),Elem(j,:)) = K(Elem(j,:),Elem(j,:)) + ST;
                    M(Elem(j,:),Elem(j,:)) = M(Elem(j,:),Elem(j,:)) + Mass;  
                    A1(Elem(j,:),Elem(j,:)) = A1(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                    [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                    A2(Elem(j,:),Elem(j,:)) = A2(Elem(j,:),Elem(j,:)) + (ysxev + ysyev)*Mass;
                    Nz(Elem(j,:),Elem(j,:)) = Nz(Elem(j,:),Elem(j,:)) + Nh;
                    Vz(Elem(j,:),Elem(j,:)) = Vz(Elem(j,:),Elem(j,:)) + Nh + V;
                end

                Awh = -eta*K(FreeNodes,FreeNodes) - A1(FreeNodes,FreeNodes) - A2(FreeNodes,FreeNodes) + ...
                      shift*M(FreeNodes,FreeNodes);
                G = M(FreeNodes,FreeNodes) - dt*Awh - dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes);
                F = -( M(FreeNodes,FreeNodes)-dt*Awh - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes))*z_temp1(FreeNodes) ...
                    + M(FreeNodes,FreeNodes)*z_temp(FreeNodes);
                zh = zeros(length(FullNodes),1);
                zh(FreeNodes) = G\F;
        
                error = Norm(Coord, Elem, zh);
                z_temp1(FreeNodes) = z_temp1(FreeNodes) + zh(FreeNodes);
            end
            Ez(k,nl) = Norm(Coord, Elem, z_temp1);
        end

        % BDF2
        for k = 3:T/dt+1
            z_iguess = z_temp1;
            error = Norm(Coord, Elem, z_temp);
            while error > 0.0001
                M = sparse(size(Coord,1),size(Coord,1)); 
                K = sparse(size(Coord,1),size(Coord,1));
                A1 = sparse(size(Coord,1),size(Coord,1));
                A2 = sparse(size(Coord,1),size(Coord,1));
                Nz = sparse(size(Coord,1),size(Coord,1));
                Vz = sparse(size(Coord,1),size(Coord,1));

                for j = 1:size(Elem,1)
                    zt = z_iguess(Elem(j,:));
                    [ST,Mass,Ah1,Nh,V] = stima_all(Coord(Elem(j,:),:),zt);
                    K(Elem(j,:),Elem(j,:)) = K(Elem(j,:),Elem(j,:)) + ST;
                    M(Elem(j,:),Elem(j,:)) = M(Elem(j,:),Elem(j,:)) + Mass;  
                    A1(Elem(j,:),Elem(j,:)) = A1(Elem(j,:),Elem(j,:)) + ys([sum(Coord(Elem(j,:),:))/3,0])*Ah1; 
                    [ysxev, ysyev] = ysxe([sum(Coord(Elem(j,:),:))/3,0]);
                    A2(Elem(j,:),Elem(j,:)) = A2(Elem(j,:),Elem(j,:)) + (ysxev + ysyev)*Mass;
                    Nz(Elem(j,:),Elem(j,:)) = Nz(Elem(j,:),Elem(j,:)) + Nh;
                    Vz(Elem(j,:),Elem(j,:)) = Vz(Elem(j,:),Elem(j,:)) + Nh + V;
                end

                Awh = -eta*K(FreeNodes,FreeNodes) - A1(FreeNodes,FreeNodes) - A2(FreeNodes,FreeNodes) + ...
                      shift*M(FreeNodes,FreeNodes);
                G = M(FreeNodes,FreeNodes) - dt*Awh - dt*exp(-shift*k*dt)*Vz(FreeNodes,FreeNodes);
                F = -( 1.5*M(FreeNodes,FreeNodes)-dt*Awh - dt*exp(-shift*k*dt)*Nz(FreeNodes,FreeNodes))*z_iguess(FreeNodes) ...
                    + 2*M(FreeNodes,FreeNodes)*z_temp1(FreeNodes) - 0.5*M(FreeNodes,FreeNodes)*z_temp(FreeNodes);
                zh = zeros(length(FullNodes),1);
                zh(FreeNodes) = G\F;
        
                error = Norm(Coord, Elem, zh);
                z_iguess(FreeNodes) = z_iguess(FreeNodes) + zh(FreeNodes);
            end
            z_temp = z_temp1;
            z_temp1 = z_iguess;
            Ez(k,nl) = Norm(Coord, Elem, z_iguess);
        end

        if (nl > 1)
            [L2ez(nl-1), H1ez(nl-1)] = Err_gaussiantest(Coord, Coord1, Elem1, ed2el, n2ed, tempz, z_iguess);
        end
        
        tempz = z_iguess;
        Coord1 = Coord;
        Elem1 = Elem; 
    end

    %% ===================== POST-PROCESSING =====================
    if abs(T - 1) < 1e-12
        Time = 0:0.001:T;
        ee = 1:1:T/dt+1;
        figure;
        plot(Time, Ez(ee,1), 'g-', 'LineWidth', 1.3); hold on
        plot(Time, Ez(ee,2), 'r-', 'LineWidth', 1.3);
        plot(Time, Ez(ee,3), 'b-', 'LineWidth', 1.3);
        xlabel('Time'); ylabel('||z||');
        legend('Level 1','Level 2','Level 3');
        title('Evolution of Ez(t) at T = 1');
        grid on;
    end

    if abs(T - 0.1) < 1e-12
        for pl = 1:(nl-2)
            oczl2(pl+1) = log(L2ez(pl)/L2ez(pl+1)) / log(h(pl)/h(pl+1));
            oczh1(pl+1) = log(H1ez(pl)/H1ez(pl+1)) / log(h(pl)/h(pl+1));
        end
        oczl2(1) = NaN;
        oczh1(1) = NaN;
        Order_of_Convergence = table(L2ez', oczl2', H1ez', oczh1');
        Order_of_Convergence.Properties.VariableNames = {'L2_Err_z','L2_Conv_z','H1_Err_z','H1_Conv_z'};
        disp('==== Order of Convergence Table (T = 0.1) ====')
        disp(Order_of_Convergence)
    end
end
