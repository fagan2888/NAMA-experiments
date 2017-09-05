close all;
clear;

rng(0, 'twister');

%% Select solvers to run

solvers = { ...
  'gpad', ...
  'nama', ...
  'gpad-scaled', ...
  'nama-scaled', ...
  'qpoases', ...
  'qpoases-ws', ...
  'mosek', ...
  'ecos', ...
  'sdpt3', ...
  'sedumi', ...
};

%% Generate problem

% Generate problem matrices (see afti16.m)

gen_afti16;

% Problem parameters

N = 50;
Q = diag([1e-4, 1e2, 1e-3, 1e2]);
R = 1e-2*eye(2);

% Define MPC problem structure

L_x = C;
L_u = speye(2);

xmin = [-0.5; -100]; xmax = [+0.5; +100];
umin = [-25; -25]; umax = [+25; +25];

mpc_prob.L_s = sparse(blkdiag(L_x, L_u));
mpc_prob.L_N = L_x;

mpc_prob.s_min = [xmin; umin];
mpc_prob.s_max = [xmax; umax];
mpc_prob.x_N_min = xmin;
mpc_prob.x_N_max = xmax;

mpc_prob.stage_w = [1e6; 1e6; inf; inf];
mpc_prob.final_w = [1e6; 1e6];

mpc_prob.Q = Q;
mpc_prob.R = R;
mpc_prob.Q_N = 100*Q;
mpc_prob.A = A;
mpc_prob.B = B;
mpc_prob.N = N;

mpc_prob.Ts = Ts;

%% Simulate system

T = 4.0; % total simulation time in seconds
t_ref = [0, 2, Inf];
x_ref = [ [0; 0; 0; 10], [0; 0; 0; 0] ];
x0 = [0; 0; 0; 0];

names = {};
times = {};
iters = {};
fops = {};
gops = {};

for k = 1:length(solvers)
    [x_sim, times_new, iters_new, fops_new, gops_new, status] = ...
        mpc_sim(mpc_prob, x0, T, solvers{k}, t_ref, x_ref);
    if status == 0
        names{end+1} = solvers{k};
        times{end+1} = times_new;
        iters{end+1} = iters_new;
        fops{end+1} = fops_new;
        gops{end+1} = gops_new;
    end
end

fprintf('%3s%12s%12s%12s%12s%12s%12s%12s%12s\n', 'id', 'avg_it', 'max_it', 'avg_f', 'max_f', 'avg_g', 'max_g', 'avg_cpu', 'max_cpu');

for k = 1:length(names)
    avg_it = mean(iters{k}(2:end));
    max_it = max(iters{k}(2:end));
    avg_cpu = mean(times{k}(2:end))*1000;
    max_cpu = max(times{k}(2:end))*1000;
    avg_fops = mean(fops{k}(2:end));
    max_fops = max(fops{k}(2:end));
    avg_gops = mean(gops{k}(2:end));
    max_gops = max(gops{k}(2:end));
    fprintf('%13s%12.2f%12d%12.2f%12d%12.2f%12d%12.3f%12.3f\n', names{k}, avg_it, max_it, avg_fops, max_fops, avg_gops, max_gops, avg_cpu, max_cpu);
end
