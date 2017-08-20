close all;
clear;

rng(0, 'twister');

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

times = {};
iters = {};

solvers = 1:10;

for id = solvers
    [x_sim, times{id}, iters{id}] = mpc_sim(mpc_prob, x0, T, id, t_ref, x_ref);
end

fprintf('%3s%12s%12s%12s%12s\n', 'id', 'avg_it', 'max_it', 'avg_cpu', 'max_cpu');

for id = solvers
    avg_it = mean(iters{id}(2:end));
    max_it = max(iters{id}(2:end));
    avg_cpu = mean(times{id}(2:end));
    max_cpu = max(times{id}(2:end));
    fprintf('%3d%12.2f%12d%12.3f%12.3f\n', id, avg_it, max_it, avg_cpu, max_cpu);
end

%% Plot results

% N_sim = size(x_sim, 2);

% figure;

% subplot(2, 1, 1);
% plot(Ts*(0:N_sim-1), x_sim(2,:));
% legend('attack angle');

% subplot(2, 1, 2);
% plot(Ts*(0:N_sim-1), x_sim(4,:));
% legend('pitch angle');

% subplot(3, 1, 3);
% semilogy(Ts*(0:N_sim-1), times{1}); hold on
% semilogy(Ts*(0:N_sim-1), times{2});
% semilogy(Ts*[0, N_sim-1], [Ts, Ts], ':');
% legend(names{1:2}, 'sampling time');

% for id = 1:length(names)
%     fprintf('%30s %8.2f %8d %7.2f %7.2f\n', names{id}, mean(iters{id}(2:end)), max(iters{id}(2:end)), mean(times{id}(2:end))*1000, max(times{id}(2:end))*1000);
% end