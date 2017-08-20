close all;
clear;

rng(0, 'twister');

ncontrols = 8;

[A, B, ts] = gen_masses(ncontrols);

n_x = size(B, 1);
n_u = size(B, 2);

% MPC parameters
Q = eye(n_x);
R = eye(n_u);

% State and input bounds
Xmax = 4; Umax = 0.5;
xmin = -Xmax*ones(n_x,1);
xmax = Xmax*ones(n_x,1);
umin = -Umax*ones(n_u,1);
umax = Umax*ones(n_u,1);

% Call LQR to get solution of Riccati equation
% and compute level set for final ellipsoidal constraint
[K, S, ~] = dlqr(A, B, Q, R);
sigma_K = svd(K);
beta = Umax/sigma_K(1);
P = chol(S); % S = P'*P
sigma_P = svd(P);
alpha = min(Xmax, beta)/sigma_P(1);

%% MPC problem structure

mpc_prob.Q = Q;
mpc_prob.R = R;
mpc_prob.Q_N = Q;
mpc_prob.A = A;
mpc_prob.B = B;
mpc_prob.L_s = sparse(blkdiag(speye(n_x), speye(n_u)));
mpc_prob.s_min = [xmin; umin];
mpc_prob.s_max = [xmax; umax];
mpc_prob.stage_w = inf(n_x + n_u, 1);
mpc_prob.x_N_ellipse = {P, alpha};

%% Set up ForBES options

opt.maxit = 100000;
opt.tol = 1e-4;
opt.display = 0;
opt.report = 0;
opt.prescale = 0;
opt.memory = 20;
opt_fama = opt; opt_fama.solver = 'fbs'; opt_fama.variant = 'fast';
opt_nama = opt; opt_nama.solver = 'nama'; opt_nama.method = 'lbfgs';

% high precision options

% opt_star = opt;
% opt_star.solver = 'nama';
% opt_star.method = 'lbfgs';
% opt_star.tol = 1e-6;

%% Solve many problems

% for k_solver = 1:5
names = {};
time_mean = {}; % average CPU time for each solver/horizon length
time_max = {}; % max CPU time for each solver/horizon length
% end
% failed = 0;

Ns = 10:10:50; % horizon lengths
nprob = 50; % number of random problems (i.e., random initial states)
nsolvers = 5;

times = zeros(length(Ns), nprob, nsolvers);
iters = zeros(length(Ns), nprob, nsolvers);

for k_N = 1:length(Ns)
    
    N = Ns(k_N);
    fprintf('N = %d\n', N);
    mpc_prob.N = N;
    
    for k_prob = 1:nprob
        
        fprintf('.');

        % Compute initial state by solving a CP so that a feasible sequence of
        % inputs/states exists.

        % restrict state/input box constraints a bit
        
        s_x = 0.4*rand(N,1); 
        s_u = 0.05*rand(N,1);
        c = randn(n_x, 1);
        
        % use ECOS to compute random initial state

        cvx_begin quiet
            cvx_solver ecos;
            cvx_precision high;
            variable x_init(n_x, 1);
            variable x_ecos(n_x, N+1);
            variable u_ecos(n_u, N);
            minimize( c'*x_init );
            subject to
                x_init == x_ecos(:,1);
                for t=1:N
                    x_ecos(:,t+1) == A*x_ecos(:,t) + B*u_ecos(:,t);
                    xmin+s_x(t) <= x_ecos(:,t) <= xmax-s_x(t);
                    umin+s_u(t) <= u_ecos(:,t) <= umax-s_u(t);
                end
                0.5*quad_form(x_ecos(:,N+1), S) <= 0.99*alpha;
        cvx_end

        % set initial state
        x0 = x_ecos(:,1);
        mpc_prob.x0 = x0;
        
        out_fama = forbes_linear_mpc(mpc_prob, opt_fama);
        names{1} = 'Fast AMA';
        t_fama = out_fama.forbes.solver.time;
        it_fama = out_fama.forbes.solver.iterations;
        times(k_N, k_prob, 1) = t_fama;
        iters(k_N, k_prob, 1) = it_fama;
        
        out_nama = forbes_linear_mpc(mpc_prob, opt_nama);
        names{2} = 'NAMA (L-BFGS)';
        t_nama = out_nama.forbes.solver.time;
        it_nama = out_nama.forbes.solver.iterations;
        times(k_N, k_prob, 2) = t_nama;
        iters(k_N, k_prob, 2) = it_nama;

        cvx_tic;
        cvx_begin quiet
            cvx_solver ecos;
            variable x_ecos(n_x,N+1);
            variable u_ecos(n_u,N);
            minimize( 0.5*(x_ecos(:)'*x_ecos(:) + u_ecos(:)'*u_ecos(:)) );
            subject to
                x_ecos(:,1) == x0;
                for t=1:N
                    x_ecos(:,t+1) == A*x_ecos(:,t) + B*u_ecos(:,t);
                    xmin <= x_ecos(:,t) <= xmax;
                    umin <= u_ecos(:,t) <= umax;
                end
                0.5*quad_form(x_ecos(:,N+1), S) <= alpha;
        cvx_end
        t_cvx = cvx_toc;
        t_ecos = t_cvx(end);
        times(k_N, k_prob, 3) = t_ecos;
        names{3} = 'ECOS (via CVX)';

        cvx_tic;
        cvx_begin quiet
            cvx_solver sdpt3;
            variable x_sdpt3(n_x,N+1);
            variable u_sdpt3(n_u,N);
            minimize( 0.5*(x_sdpt3(:)'*x_sdpt3(:) + u_sdpt3(:)'*u_sdpt3(:)) );
            subject to
                x_sdpt3(:,1) == x0;
                for t=1:N
                    x_sdpt3(:,t+1) == A*x_sdpt3(:,t) + B*u_sdpt3(:,t);
                    xmin <= x_sdpt3(:,t) <= xmax;
                    umin <= u_sdpt3(:,t) <= umax;
                end
                0.5*quad_form(x_sdpt3(:,N+1), S) <= alpha;
        cvx_end
        t_cvx = cvx_toc;
        t_sdpt3 = t_cvx(end);
        times(k_N, k_prob, 4) = t_sdpt3;
        names{4} = 'SDPT3 (via CVX)';

        cvx_tic;
        cvx_begin quiet
            cvx_solver sedumi;
            variable x_sedumi(n_x,N+1);
            variable u_sedumi(n_u,N);
            minimize( 0.5*(x_sedumi(:)'*x_sedumi(:) + u_sedumi(:)'*u_sedumi(:)) );
            subject to
                x_sedumi(:,1) == x0;
                for t=1:N
                    x_sedumi(:,t+1) == A*x_sedumi(:,t) + B*u_sedumi(:,t);
                    xmin <= x_sedumi(:,t) <= xmax;
                    umin <= u_sedumi(:,t) <= umax;
                end
                0.5*quad_form(x_sedumi(:,N+1), S) <= alpha;
        cvx_end
        t_cvx = cvx_toc;
        t_sedumi = t_cvx(end);
        times(k_N, k_prob, 5) = t_sedumi;
        names{5} = 'SeDuMi (via CVX)';

    end

    for k_solver = 1:length(names)
        fprintf('\n%25s avg %7.2f max %7.2f', names{k_solver}, mean(times(k_N, :, k_solver)), max(times(k_N, :, k_solver)));
    end
    fprintf('\n');

end

%% Plot results

subplot(1, 2, 1);

for k_solver = 1:length(names)
    semilogy(Ns, mean(times(:, :, k_solver), 2), 'LineWidth', 2); hold on;
end
legend(names{:});
str_title = ['Average CPU time (s)'];
title(str_title);
xlabel('Control horizon');
ylabel('Average CPU time (s)');
grid on;

subplot(1, 2, 2);

for k_solver = 1:length(names)
    semilogy(Ns, max(times(:, :, k_solver), [], 2), 'LineWidth', 2); hold on;
end
legend(names{:});
str_title = ['Max CPU time (s)'];
title(str_title);
xlabel('Control horizon');
ylabel('Maximum CPU time (s)');
grid on;
