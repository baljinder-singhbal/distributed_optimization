% In the proper folder there's the random_seed_task3 used for the experiment, %
% open it in the workspace to repeat the experiment, otherwise use the  %
% variable keep_prev_exp to control the generation of initial data.     %

clearvars -except random_seed_task3;
close all

%% CUSTOM PARAMETERS %%

    % Probability related to the adjacency matrix
    adj_p = 0.3;
    
    % Nodes numbers
    N_gen = 4;
    N_stor = 3;
    N_conl = 2;
    N_tr = 1;
    
    % Period
    T = 11;
    
    % Max iterations of the algorithm
    MAXITERS = 1500;

    % Step size
    eta=0.002;

    % Error Threshold
    threshold = 0.000001;

    % Repeatability (1 = repeat, 0 = new experiment)
    keep_prev_exp = 1;
    
    % Repeatability function
        if keep_prev_exp~=0 && exist('random_seed_task3','var')
            rng(random_seed_task3);
        else
            random_seed_task3 = rng;
        end

    
%% PROBLEM DEFINITION %%


    %% Generator Nodes %%
    
        % Constraints values
        p_gen_sup = 100*rand(T+1,1);
        p_gen_inf = -100*rand(T+1,1);
        B = (-eye(T+1) + diag(ones(T+1-abs(1),1),1));
        B(end,:) = [];
        r_gen_sup = 100*rand(T,1);
        r_gen_inf = -100*rand(T,1);
        
        % Constraints composition for fmincon
        gen_ineq_matrix = [-B; B];
        gen_ineq_bounds = [-r_gen_inf; r_gen_sup];
        
        % Node Function
        alpha_1 = rand(1);
        alpha_2 = rand(1);
        f_gen = @(p_gen)alpha_2 * p_gen'*p_gen + alpha_1* p_gen'*(ones(T+1,1));
        
    %% Storage Nodes %%
        
        % Constraints values
        c_stor = 100*rand(T+1,1);
        d_stor = 100*rand(T+1,1);
        Tr = tril(ones(T,T+1));
        q_stor_max = 500*rand(T,1);
        
        % Constraints composition for fmincon
        stor_ineq_matrix = [-Tr; Tr];
        stor_ineq_bounds = [zeros(T,1); q_stor_max];
        
        % Node Function
        eps = rand(1);
        f_stor = @(p_stor)eps*(p_stor')*p_stor;

    %% Controllable Load Nodes %%
    
        % Constraints values
        P_conl = 100*rand(T+1,1);
        p_des=rand(T+1,1);
        
        % Node Function
        beta = rand(1);
        f_conl = @(p_conl)beta*max(zeros(T+1,1),p_des-p_conl)'*ones(T+1,1) + eps*p_conl'*p_conl;
        
    %% Trade Nodes %%
    
        % Constraints values
        E_tr = 100*rand(T+1,1);
        D = rand(T+1,1);
        
        % Node Function
        c_1 = rand(1);
        c_2 = rand(1);
        f_tr = @(p_tr)(-c_1*p_tr+c_2*abs(p_tr))'*ones(T+1,1)+eps*(p_tr')*p_tr;

        
%% VARIABLES AND FUNCTIONS %%

    %% Nodes %%
        
        % Number of agents
        N = N_gen + N_stor + N_conl + N_tr;
        
        % Indexex Calculation
        i_gen = 1;
        i_stor = i_gen + N_gen;
        i_conl = i_stor + N_stor;
        i_tr = i_conl + N_conl;
    
    %% Map %%

        % Adjacency matrix generation
        while 1
            
          Adj = binornd(1,adj_p,N,N);
          I_NN = eye(N,N);
          notI = ~I_NN;
          Adj = Adj.*notI;
          Adj = or(Adj,Adj');
          test = (I_NN+Adj)^N;
          
          if ~any(any(~test))
            fprintf('the graph is connected\n');
            break
          else
            fprintf('the graph is not connected\n');
          end
          
        end
        
        % Weighted consensus matrix generation
        WW = zeros(N,N);
        DEGREE = sum(Adj); %Degree of the WW matrix
        for i = 1:N
          N_i = find(Adj(:,i) == 1)';
          for j = N_i
            WW(i,j) = 1/(1 + max(DEGREE(i),DEGREE(j) ));
          end
          WW(i,i) = 1 - sum(WW(i,:));
        end

    %% Algorithm %%
        
        % f_i matrices
        f = zeros(1,MAXITERS);

        % P_error
        P_error = zeros(T+1,MAXITERS);

        % P matrices
        P = zeros(T+1,MAXITERS,N);
        
        % Lambda matrices
        LM = zeros(T+1,MAXITERS,N);
        %LM(:,1,:) = -rand(T+1,N);

        % Gradient matrices
        grad_q = zeros(T+1,MAXITERS,N);

        % Gradient descent matrices
        S = zeros(T+1,MAXITERS,N);
        
        
     %% Other options %%
     
        % fmincon options
        options = optimoptions('fmincon','Display','off', 'Algorithm', 'sqp'); 

%% CENTRALIZED FUNCTION GENERATION %%
    f_c = @(p) 0;

        for i = 1:N

                    if (i_gen <= i) && (i <= (i_stor-1))
                        f_c = @(p) f_c(p) + f_gen(p(:,i));
                    end

                    if (i_stor <= i) && (i <= (i_conl-1))
                        f_c = @(p) f_c(p) + f_stor(p(:,i));
                    end

                    if (i_conl <= i) && (i <= (i_tr-1))
                        f_c = @(p) f_c(p) + f_conl(p(:,i));
                    end

                    if (i_tr <= i) && (i <= N)
                        f_c = @(p) f_c(p) + f_tr(p(:,i));
                    end

        end
        
%% MULTIOBJECTIVE CENTRALIZED OPTIMIZATION %%


    cvx_begin quiet
       variables p(T+1, N)
       dual variable lambda
       minimize(f_c(p))
       subject to
            for i = 1:N
    
                if (i_gen <= i) && (i <= (i_stor-1))
                    p(:, i) <= p_gen_sup
                    p(:, i) >= p_gen_inf
                    B*p(:, i) >= r_gen_inf
                    B*p(:, i) <= r_gen_sup
                end

                if (i_stor <= i) && (i <= (i_conl-1))
                    p(:, i) <= c_stor
                    p(:, i) >= -d_stor
                    Tr*p(:, i) >= 0
                    Tr*p(:, i) <= q_stor_max
                end

                if (i_conl <= i) && (i <= (i_tr-1))
                    p(:, i) <= P_conl
                    p(:, i) >= -P_conl
                end

                if (i_tr <= i) && (i <= N)
                    p(:, i) <=  E_tr
                    p(:, i) >= -E_tr
                end

            end

            lambda : D == sum(p')'
    cvx_end
    
  p_opt = p;
  
  f_star = f_c(p_opt);
  
  
%% INITIALIZATION OF L, P AND S

            
        for i = 1:N            

            if (i_gen <= i) && (i <= (i_stor-1))
                P(:,1,i) = fmincon(@(p_gen) f_gen(p_gen) + LM(:,1,i)'*p_gen,LM(:,1,i),gen_ineq_matrix,gen_ineq_bounds,[],[],p_gen_inf,p_gen_sup,'',options);
                S(:,1,i)= -N*P(:,1,i);
            end

            if (i_stor <= i) && (i <= (i_conl-1))
                P(:,1,i) = fmincon(@(p_stor) f_stor(p_stor) + LM(:,1,i)'*p_stor,LM(:,1,i),stor_ineq_matrix,stor_ineq_bounds,[],[],-d_stor,c_stor,'',options);
                S(:,1,i)= -N*P(:,1,i);
            end

            if (i_conl <= i) && (i <= (i_tr-1))
                P(:,1,i) = fmincon(@(p_conl) f_conl(p_conl) + LM(:,1,i)'*p_conl,LM(:,1,i),[],[],[],[],-P_conl,P_conl,'',options);
                S(:,1,i)= -N*P(:,1,i);
            end

            if (i_tr <= i) && (i <= N)
                P(:,1,i) = fmincon(@(p_tr) f_tr(p_tr) + LM(:,1,i)'*(p_tr-D),LM(:,1,i),[],[],[],[],-E_tr,E_tr,'',options);
                S(:,1,i)= -N*(P(:,1,i)-D);
            end
           
        end
        

%% MULTIOBJECTIVE DISTRIBUTED OPTIMIZATION %%
    

    % Waitbar initialization
    message = sprintf('f - f* is: \n unknown  \n max_i ||λ^t_i - λ^*|| is: \n unknown \n');
    w = waitbar(0,message,'Name','Optimizing...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    w_c = get(w, 'children');
    set(w_c(1), 'String', 'Stop');
    setappdata(w,'canceling',0);
    
    % Time iterations
    for t = 1:MAXITERS-1
                
        % Executing the algorithm for each agent
        for i = 1:N
            
            % Update of lambdas wrt to the neighbours
            for j = 1:N
              LM(:,t+1,i) = LM(:,t+1,i) + WW(i,j)*LM(:,t,j);
              S(:,t+1,i) =  S(:,t+1,i) + WW(i,j)*S(:,t,j);
            end

            % Gradient descent application
            LM(:,t+1,i) = LM(:,t+1,i) - eta*S(:,t,i);
                           
            % 1. Local minimization
            % 2. Calculation of the full function value
            % 3. Calculation of the p error
            
            if (i_gen <= i) && (i <= (i_stor-1))
                P(:,t+1,i) = fmincon(@(p_gen) f_gen(p_gen) + LM(:,t,i)'*(p_gen),P(:,1,i),gen_ineq_matrix,gen_ineq_bounds,[],[],p_gen_inf,p_gen_sup,'',options);
                f(1,t+1) = f(:,t+1) + f_gen(P(:,t+1,i));
                P_error(:,t,i) = p_opt(:,i) - P(:,t,i);
            end

            if (i_stor <= i) && (i <= (i_conl-1))
                P(:,t+1,i) = fmincon(@(p_stor) f_stor(p_stor) + LM(:,t,i)'*(p_stor),P(:,1,i),stor_ineq_matrix,stor_ineq_bounds,[],[],-d_stor,c_stor,'',options);
                f(1,t+1) = f(:,t+1) + f_stor(P(:,t+1,i));
                P_error(:,t,i) = p_opt(:,i) - P(:,t,i);
            end

            if (i_conl <= i) && (i <= (i_tr-1))
                P(:,t+1,i) = fmincon(@(p_conl) f_conl(p_conl) + LM(:,t,i)'*(p_conl),P(:,1,i),[],[],[],[],-P_conl,P_conl,'',options);
                f(1,t+1) = f(:,t+1) + f_conl(P(:,t+1,i));
                P_error(:,t,i) = p_opt(:,i) - P(:,t,i);
            end

            if (i_tr <= i) && (i <= N)
                P(:,t+1,i) = fmincon(@(p_tr) f_tr(p_tr) + LM(:,t,i)'*(p_tr-D),P(:,1,i),[],[],[],[],-E_tr,E_tr,'',options);
                f(1,t+1) = f(:,t+1) + f_tr(P(:,t+1,i));
                P_error(:,t,i) = p_opt(:,i) - P(:,t,i);
            end
            
            
            % Gradients computation
            grad_q(:, t,i) = -N*P(:,t,i);
            grad_q(:, t+1,i) = -N*P(:,t+1,i);
            
            % Gradient descent update
            S(:,t+1,i) = S(:,t+1,i) + (grad_q(:,t+1,i) - grad_q(:,t,i));

        end
        
        % Current f-f* error
        f_inst_error = f(1,t+1) - f_star;
        
        % Current max_i ||λ^t_i - λ^*||
        lm_inst_error = max(vecnorm(LM(:,t,:) - lambda),[],3);
        
        message = sprintf('f - f* is: \n %.7f  \n max_i ||λ^t_i - λ^*|| is: \n %.7f \n',f_inst_error,lm_inst_error);
        
        % Show waitbar to monitor the progress
        waitbar(t/MAXITERS,w,message);

        
        % Exit if the error is lower than a threshold or cancel is clicked
        if (t>=5 & abs(lm_inst_error) <= threshold) | (getappdata(w,'canceling'))
           MAXITERS = t;
           break;
        end
        
    end
    
    % Delete waitbar object
    delete(w)
    
        
%% PLOTS %%

    % Create the instances of the group and the figures
    
        desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
        myGroup = desktop.addGroup('Plots');
        desktop.setGroupDocked('Plots', 0);
        myDim   = java.awt.Dimension(5, 1);
        desktop.setDocumentArrangement('Plots', 1, myDim);
        figures = gobjects(1, 5);
        bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        set(0,'defaultAxesFontSize',14);
        
    % | f - f* |
        
        % Create the figure
        figures(1) = figure('WindowStyle', 'docked','Name', sprintf('| f - f* |'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(1)), 'javaframe'), 'GroupName', 'Plots');
        
        % Compute the error | f - f* |
        f_error_plot = f(1,(1:MAXITERS-2))-f_star*ones(1,MAXITERS-2);
        
        % Plot
        semilogy(1:(MAXITERS-2), abs(f_error_plot),'b','LineWidth',1.5);
        xlabel('t');
        grid on
        
        legend('   | f - f* |','FontSize',14);
        legend('boxoff');

        title('| f - f* |');
        
     % max_i ||λ^t_i - λ^*||
    
        % Create the figure
        figures(2) = figure('WindowStyle', 'docked','Name', sprintf('max_i ||λ^t_i - λ^*||'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(2)), 'javaframe'), 'GroupName', 'Plots');
        
        % Compute max_i ||λ^t_i - λ^*||
        lm_error_plot = max(vecnorm(LM(:,1:(MAXITERS-2),:) - lambda*ones(1,MAXITERS-2)),[],3);
        
        % Plot
        semilogy(1:(MAXITERS-2), lm_error_plot,'r','LineWidth',1.5);
        xlabel('t') ;
        grid on
        
        legend('   max_i ||λ^t_i - λ^*||','FontSize',14);
        legend('boxoff');
        
        title('max_i ||λ^t_i - λ^*||');
        
    % max_i || S^t_i ||
    
        % Create the figure
        figures(3) = figure('WindowStyle', 'docked','Name', sprintf('max_i || S^t_i ||'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(3)), 'javaframe'), 'GroupName', 'Plots');
        
        % Compute max_i || S^t_i ||
        s_vector_plot = max(vecnorm(S(:,1:(MAXITERS-2),:)),[],3);
        
        % Plot
        semilogy(1:MAXITERS-2, s_vector_plot,'g','LineWidth',1.5);
        xlabel('t');
        grid on
        legend('   max_i || S^t_i ||','FontSize',14)
        legend('boxoff');
        
        title('max_i || S^t_i ||');
        
	% max_i || P^t_i - P* ||
    
        % Create the figure
        figures(4) = figure('WindowStyle', 'docked','Name', sprintf('max_i || P^t_i - P* ||'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(4)), 'javaframe'), 'GroupName', 'Plots');
        
        % Compute max_i || P^t_i - P* ||
        p_error_plot = max(vecnorm(P_error(:,1:MAXITERS-2,:)),[],3);
        
        % Plot
        semilogy(1:MAXITERS-2, p_error_plot,'k','LineWidth',1.5);
        xlabel('t') ;
        grid on
        legend('   max_i || P^t_i - P* ||','FontSize',14);
        legend('boxoff');
        
        title('max_i || P^t_i - P* ||');
        
     % max_i ||λ^t_i - λ_{mean}||
    
        % Create the figure
        figures(5) = figure('WindowStyle', 'docked','Name', sprintf('max_i ||λ^t_i - λ_{mean}||'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(5)), 'javaframe'), 'GroupName', 'Plots');
        
        % Compute max_i ||λ^t_i - λ_{mean}||
        lm_error_mean_plot = max(vecnorm(LM(:,1:(MAXITERS-2),:) - mean(LM(:,1:(MAXITERS-2),:),3)),[],3);
        
        % Plot
        semilogy(1:(MAXITERS-2), lm_error_mean_plot,'m','LineWidth',1.5);
        xlabel('t');
        grid on
        
        legend(' max_i ||λ^t_i - λ_{mean}||','FontSize',14);
        legend('boxoff');
        
        title('max_i ||λ^t_i - λ_{mean}||');
        hold on
        
    % Node Map
    
        % Create the figure
        figures(6) = figure('WindowStyle', 'docked','Name', sprintf('Node Map'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(6)), 'javaframe'), 'GroupName', 'Plots');
        
        % Plot the map
        G = graph(WW,'omitselfloops');
        LWidths = G.Edges.Weight/max(G.Edges.Weight);
        h = plot(G,'LineWidth',LWidths);
        title('Node Map');
        
    % Make the window appears
    drawnow;
    
    % Remove warning of obsolete javaframe
    warning(bakWarn);