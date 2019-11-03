% In the proper folder there's the random_seed_task1 used for the experiment, %
% open it in the workspace to repeat the experiment, otherwise use the  %
% variable keep_prev_exp to control the generation of initial data.     %

clearvars -except random_seed_task1;
close all

%% CUSTOM PARAMETERS %%

    % Probability related to the adjacency matrix
    adj_p = 0.3;
    
    % Nodes numbers
    N = 5;
    
    % Dimension of the state of each agent
    T = 3;
    
    % Max iterations of the algorithm
    MAXITERS = 10000;

    % Step size
    eta=0.001;

    % Error Threshold
    threshold = 0.0000000000001;

    % Repeatability (1 = repeat, 0 = new experiment)
    keep_prev_exp = 1;
    
    % Repeatability function
        if keep_prev_exp~=0 && exist('random_seed_task1','var')
            rng(random_seed_task1);
        else
            random_seed_task1 = rng;
        end

    
%% PROBLEM DEFINITION %%

    % Matrices for the function %
        
        Q = bsxfun(@times, eye(T), rand(T,T,N)); % Diagonal matrices
        R = rand(T,N);
        
    % Map %

        % Adjacency matrix generation
        while 1
            
          Adj = binornd(1,adj_p,N,N);
          I_NN = eye(N,N);
          notI = ~I_NN;
          Adj = Adj.*notI; % removes element on the diagonal
          Adj = or(Adj,Adj'); % makes it symmetric
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

    % Algorithm %
        
        % q vector
        q = zeros(1,MAXITERS);
        
        % P_error
        P_error = zeros(T,MAXITERS);
        
        % Lambda matrices
        LM = zeros(T,MAXITERS,N);

        % Gradient matrices
        grad_q = zeros(T,MAXITERS,N);

        % Gradient descent matrices
        S = zeros(T,MAXITERS,N);
        

  %% OPTIMAL SOLUTION IN A NON DISTRIBUTED WAY %%

    sum_Q=zeros(T,T);
    sum_R=zeros(T,1);
    for i = 1:N
         sum_Q = sum_Q + Q(:,:,i);
         sum_R = sum_R + R(:,i);
    end

    lambda = quadprog(2*N*sum_Q,N*sum_R);
    q_star = - ( lambda'*sum_Q*lambda + sum_R'*lambda );
  
  
%% INITIALIZATION OF S

    for i = 1:N
        S(:,1,i)=N*(2*Q(:,:,i)*LM(:,1,i) + R(:,i));
    end


%% MULTIOBJECTIVE DISTRIBUTED OPTIMIZATION %%
    

    % Waitbar initialization
    message = sprintf('q - q* is: \n unknown  \n max_i ||λ^t_i - λ^*|| is: \n unknown \n');
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
            
            % Gradients computation
            grad_q(:, t,i) = (2*Q(:,:,i)*LM(:,t,i) + R(:,i));
            grad_q(:, t+1,i) = (2*Q(:,:,i)*LM(:,t+1,i) + R(:,i));
            
            % Gradient descent update
            S(:,t+1,i) = S(:,t+1,i) + N*(grad_q(:,t+1,i) - grad_q(:,t,i));
            
            q(1,t+1) = q(1,t+1) - ( LM(:,t+1,i)'*Q(:,:,i)*LM(:,t+1,i) + R(:,i)'*LM(:,t+1,i) );

        end
        
        
        
        % Current q-q* error
        q_inst_error = q(1,t+1) - q_star;
        
        % Current max_i ||λ^t_i - λ^*||
        lm_inst_error = max(vecnorm(LM(:,t,:) - lambda),[],3);
        
        message = sprintf('q - q* is: \n %.7f  \n max_i ||λ^t_i - λ^*|| is: \n %.7f \n',q_inst_error,lm_inst_error);
        
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
        
    % | q - q* |
        
        % Create the figure
        figures(1) = figure('WindowStyle', 'docked','Name', sprintf('| q - q* |'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(1)), 'javaframe'), 'GroupName', 'Plots');
        
        % Compute the error | q - q* |
        q_error_plot = q(1,(1:MAXITERS-2))-q_star*ones(1,MAXITERS-2);
        
        % Plot
        semilogy(1:(MAXITERS-2), abs(q_error_plot),'b','LineWidth',1.5);
        xlabel('t');
        grid on
        
        legend('   | q - q* |','FontSize',14);
        legend('boxoff');

        title('| q - q* |');
        
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
        
        
     % max_i ||λ^t_i - λ_{mean}||
    
        % Create the figure
        figures(4) = figure('WindowStyle', 'docked','Name', sprintf('max_i ||λ^t_i - λ_{mean}||'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(4)), 'javaframe'), 'GroupName', 'Plots');
        
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
        figures(5) = figure('WindowStyle', 'docked','Name', sprintf('Node Map'), 'NumberTitle', 'off');
        
        % Put the created figure in the group
        set(get(handle(figures(5)), 'javaframe'), 'GroupName', 'Plots');
        
        % Plot the map
        G = graph(WW);
        LWidths = G.Edges.Weight/max(G.Edges.Weight);
        h = plot(G,'LineWidth',LWidths,'EdgeLabel',G.Edges.Weight);
        title('Node Map');
        
    % Make the window appears
    drawnow;
    
    % Remove warning of obsolete javaframe
    warning(bakWarn);