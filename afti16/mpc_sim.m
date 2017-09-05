function [x_sim, times, iters, fops, gops, status] = mpc_sim(mpc_prob, x_curr, T, solver, t_ref, x_ref)
% MPC_SIM simulates the discrete-time linear dynamics whose description is
% in the mpc_prob structure, for a simulation time T, using the specified
% solver.
%
% The solver is selected as follows:
%
% solver = 'gpad'          =>   GPAD (no scaling)
%        = 'nama'          =>   NAMA (no scaling, L-BFGS)
%        = 'gpad-scaled'   =>   GPAD (Jacobi scaling)
%        = 'nama-scaled'   =>   NAMA (Jacobi scaling, L-BFGS)
%        = 'qpoases'       =>   qpOASES
%        = 'qpoases-ws'    =>   qpOASES (warm-started)
%        = 'mosek'         =>   MOSEK
%        = 'ecos'          =>   ECOS
%        = 'sdpt3'         =>   SDPT3
%        = 'sedumi'        =>   SeDuMi

N = mpc_prob.N;
Ts = mpc_prob.Ts;
N_sim = floor(T/Ts); % simulation steps
n_x = size(mpc_prob.Q, 2);
n_u = size(mpc_prob.R, 2);

times = zeros(1, N_sim); % solver time at each simulation step
iters = zeros(1, N_sim); % solver iterations at each simulation step
fops = zeros(1, N_sim); % solver operations wrt f at each simulation step
gops = zeros(1, N_sim); % solver operations wrt g at each simulation step

x_sim = zeros(n_x, 1);

k_ref = 1;

out = [];
status = 0;

fprintf('Running %13s: ', solver);

for k_sim = 1:N_sim

    % update reference if necessary
    if nargin > 4 && k_sim*Ts >= t_ref(k_ref)
        mpc_prob.xref = x_ref(:, k_ref);
        k_ref = k_ref+1;
    end

    % update initial state
    x_sim(:,k_sim) = x_curr;
    mpc_prob.x0 = x_curr;

    [H, q, A_eq, b_eq, A_ineq, l_ineq, u_ineq, A_slack, q_slack] = mpc_to_qp(mpc_prob);
    nvar = (N+1)*n_x+N*n_u;
    nslack = length(q_slack);

    % We will use the following matrices and vectors for QP solvers

    H_qp = blkdiag(H, sparse(nslack, nslack));
    q_qp = [q; q_slack];
    A_eq_qp = [A_eq, sparse(size(A_eq, 1), nslack)];
    A_ineq_qp = [A_ineq, A_slack];
    lx_qp = [-inf(nvar, 1); zeros(nslack, 1)];

    opt.maxit = 500000;
    opt.tol = 1e-4;
    opt.display = 0;
    opt.report = 0;
    opt.memory = 20;

    switch solver

    case 'gpad'

        if exist('forbes_linear_mpc') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end

        opt.solver = 'fbs';
        opt.variant = 'fast';
        opt.prescale = 0;
        out = forbes_linear_mpc(mpc_prob, opt, out);
        times(k_sim) = out.forbes.solver.time;
        iters(k_sim) = out.forbes.solver.iterations;
        fops(k_sim) = out.forbes.solver.operations.gradf1;
        gops(k_sim) = out.forbes.solver.operations.proxg;
        u_curr = out.u(:, 1);

    case 'nama'

        if exist('forbes_linear_mpc') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end

        opt.solver = 'nama';
        opt.method = 'lbfgs';
        opt.prescale = 0;
        out = forbes_linear_mpc(mpc_prob, opt, out);
        times(k_sim) = out.forbes.solver.time;
        iters(k_sim) = out.forbes.solver.iterations;
        fops(k_sim) = out.forbes.solver.operations.gradf1;
        gops(k_sim) = out.forbes.solver.operations.proxg;
        u_curr = out.u(:, 1);

    case 'gpad-scaled'

        if exist('forbes_linear_mpc') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end

        opt.solver = 'fbs';
        opt.variant = 'fast';
        opt.prescale = 1;
        out = forbes_linear_mpc(mpc_prob, opt, out);
        times(k_sim) = out.forbes.solver.time;
        iters(k_sim) = out.forbes.solver.iterations;
        fops(k_sim) = out.forbes.solver.operations.gradf1;
        gops(k_sim) = out.forbes.solver.operations.proxg;
        u_curr = out.u(:, 1);

    case 'nama-scaled'

        if exist('forbes_linear_mpc') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end

        opt.solver = 'nama';
        opt.method = 'lbfgs';
        opt.prescale = 1;
        out = forbes_linear_mpc(mpc_prob, opt, out);
        times(k_sim) = out.forbes.solver.time;
        iters(k_sim) = out.forbes.solver.iterations;
        fops(k_sim) = out.forbes.solver.operations.gradf1;
        gops(k_sim) = out.forbes.solver.operations.proxg;
        u_curr = out.u(:, 1);

    case 'qpoases'

        if exist('qpOASES') ~= 3
            fprintf('not installed');
            status = 1;
            break
        end

        t0 = tic();
        [xus_oases,fval_oases,flag_oases,iter_oases,lambda_oases] = ...
            qpOASES(H_qp, q_qp, [A_eq_qp; A_ineq_qp], lx_qp, [], [b_eq; l_ineq], [b_eq; u_ineq]);
        t_qpoases = toc(t0);
        times(k_sim) = t_qpoases;
        iters(k_sim) = 0;
%       xs(:,k_sim+1) = xus_oases(n);
        u_curr = xus_oases(n_x+1:n_x+n_u);

    case 'qpoases-ws'

        if exist('qpOASES') ~= 3
            fprintf('not installed');
            status = 1;
            break;
        end

        if k_sim == 1
            % first (cold) run
            t0 = tic();
            [QP_oasesws,xus_oasesws,fval_oasesws,flag_oasesws,iter_oasesws,lambda_oasesws] = ...
                qpOASES_sequence('i', H_qp, q_qp, [A_eq_qp; A_ineq_qp], lx_qp, [], [b_eq; l_ineq], [b_eq; u_ineq]);
            t_oasesws = toc(t0);
        else
            % warm started runs
            t0 = tic();
            [xus_oasesws,fval_oasesws,flag_oasesws,iter_oasesws,lambda_oasesws] = ...
                qpOASES_sequence('h', QP_oasesws, q_qp, lx_qp, [], [b_eq; l_ineq], [b_eq; u_ineq]);
            t_oasesws = toc(t0);
        end
        times(k_sim) = t_oasesws;
        iters(k_sim) = 0;
        u_curr = xus_oasesws(n_x+1:n_x+n_u);

    case 'mosek'

        if exist('mskqpopt') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end

        param_msk.MSK_IPAR_LOG = 0;
        t0 = tic();
        out_msk = mskqpopt(H_qp, q_qp, [A_eq_qp; A_ineq_qp], [b_eq; l_ineq], [b_eq; u_ineq], lx_qp, [], param_msk);
        t_msk = toc(t0);
        xus_msk = out_msk.sol.itr.xx(1:nvar);
        times(k_sim) = t_msk;
        iters(k_sim) = 0;
        u_curr = xus_msk(n_x+1:n_x+n_u);

    case 'ecos'

        if exist('cvx_begin') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end
        
        solver_backup = cvx_solver;
        try
            cvx_solver('ecos');
        catch
            fprintf('not installed');
            status = 1;
            break;
        end
        cvx_solver(solver_backup);

        cvx_tic;
        cvx_begin quiet
            variable xu_ecos(nvar);
            variable s_ecos(nslack);
            minimize( 0.5*(xu_ecos'*H*xu_ecos) + q'*xu_ecos + q_slack'*s_ecos);
            subject to
                A_eq*xu_ecos == b_eq;
                l_ineq <= A_ineq*xu_ecos + A_slack*s_ecos <= u_ineq;
                s_ecos >= 0;
        cvx_end
        t_cvx = cvx_toc;
        t_ecos = t_cvx(end);
        times(k_sim) = t_ecos;
        iters(k_sim) = 0;
        u_curr = xu_ecos(n_x+1:n_x+n_u);

    case 'sdpt3'

        if exist('cvx_begin') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end
        
        solver_backup = cvx_solver;
        try
            cvx_solver('sdpt3');
        catch
            fprintf('not installed');
            status = 1;
            break;
        end
        cvx_solver(solver_backup);

        cvx_tic;
        cvx_begin quiet
            variable xu_sdpt3(nvar);
            variable s_sdpt3(nslack);
            minimize( 0.5*(xu_sdpt3'*H*xu_sdpt3) + q'*xu_sdpt3 + q_slack'*s_sdpt3);
            subject to
                A_eq*xu_sdpt3 == b_eq;
                l_ineq <= A_ineq*xu_sdpt3 + A_slack*s_sdpt3 <= u_ineq;
                s_sdpt3 >= 0;
        cvx_end
        t_cvx = cvx_toc;
        t_sdpt3 = t_cvx(end);
        times(k_sim) = t_sdpt3;
        iters(k_sim) = 0;
        u_curr = xu_sdpt3(n_x+1:n_x+n_u);

    case 'sedumi'

        if exist('cvx_begin') ~= 2
            fprintf('not installed');
            status = 1;
            break;
        end
        
        solver_backup = cvx_solver;
        try
            cvx_solver('sedumi');
        catch
            fprintf('not installed');
            status = 1;
            break;
        end
        cvx_solver(solver_backup);

        cvx_tic;
        cvx_begin quiet
            variable xu_sedumi(nvar);
            variable s_sedumi(nslack);
            minimize( 0.5*(xu_sedumi'*H*xu_sedumi) + q'*xu_sedumi + q_slack'*s_sedumi);
            subject to
                A_eq*xu_sedumi == b_eq;
                l_ineq <= A_ineq*xu_sedumi + A_slack*s_sedumi <= u_ineq;
                s_sedumi >= 0;
        cvx_end
        t_cvx = cvx_toc;
        t_sedumi = t_cvx(end);
        times(k_sim) = t_sedumi;
        iters(k_sim) = 0;
        u_curr = xu_sedumi(n_x+1:n_x+n_u);

    otherwise

        fprintf('unknown solver');
        status = 1;
        break;

    end

    % evolve system
    x_curr = mpc_prob.A*x_curr + mpc_prob.B*u_curr;

    fprintf('.');

end

fprintf('\n');
