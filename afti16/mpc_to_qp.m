function [H, q, A_eq, b_eq, A_ineq, b_ineq_low, b_ineq_upp, A_slack, q_slack] = mpc_to_qp(mpc_prob)
% MPC_TO_QP transform a discrete-time linear MPC problem, whose description
% is contained in the mpc_prob structure, into a QP in the following form
% 
%   minimize    (1/2)*(x'*H*x) + q'*x + q_slack'*s
%   subject to  A_eq*x = b_eq
%               b_ineq_low <= A_ineq*x + A_slack*s <= b_ineq_upp
% 
% where the decision variables are x and s. Variable x contains states and
% inputs of the original MPC problem, while s is a slack vector and is used
% to encode soft inequality constraints (linearly penalized).

x0 = mpc_prob.x0;
N = mpc_prob.N;
Q = mpc_prob.Q;
R = mpc_prob.R;
Q_N = mpc_prob.Q_N;
A = mpc_prob.A;
B = mpc_prob.B;
L_s = mpc_prob.L_s;
smin = mpc_prob.s_min;
smax = mpc_prob.s_max;

if ~isfield(mpc_prob, 'xref')
    flag_N = 0;
else
    flag_N = 1;
    L_N = mpc_prob.L_N;
    xmin_N = mpc_prob.x_N_min;
    xmax_N = mpc_prob.x_N_max;
end

if ~isfield(mpc_prob, 'xref')
    xref = [];
else
    xref = mpc_prob.xref;
end

n_x = size(Q, 2);
n_u = size(R, 2);

n_var = (N+1)*n_x + N*n_u;

if isempty(xref), xref = zeros(n_x, 1); end

block_eq = [A, B, -eye(n_x)];
A_dyn = sparse(N*n_x, (N+1)*n_x + N*n_u);
diag_A_ineq = {}; b_ineq_upp = []; b_ineq_low = []; q = []; w = [];
diag_H = {};
for i = 0:N-1 % build up constraints
    basei = i*n_x;
    basej = i*(n_x+n_u);
    A_dyn(basei+1:basei+n_x, basej+1:basej+2*n_x+n_u) = block_eq;
    diag_A_ineq{i+1} = [L_s];
    b_ineq_upp = [b_ineq_upp; smax];
    b_ineq_low = [b_ineq_low; smin];
    diag_H{2*i+1} = Q;
    diag_H{2*i+2} = R;
    q = [q; -Q*xref; zeros(n_u, 1)];
    w = [w; mpc_prob.stage_w];
end
if flag_N
    diag_A_ineq{N+1} = [L_N];
    b_ineq_upp = [b_ineq_upp; xmax_N];
    b_ineq_low = [b_ineq_low; xmin_N];
    w = [w; mpc_prob.final_w];
end
q = [q; -Q_N*xref];
diag_H{2*N+1} = Q_N;

% add slack variables for soft constraints
block_slack = sparse([1, -1]);
flag_slack = w < inf;
mask_slack = speye(length(w));
mask_slack = mask_slack(:, flag_slack);
A_slack = kron(mask_slack, block_slack);
n_s = size(A_slack, 2);
q_slack = reshape([w(flag_slack), w(flag_slack)]', n_s, 1);

H = sparse(blkdiag(diag_H{:}));
A_ineq = sparse(blkdiag(diag_A_ineq{:}));
if flag_N == 0
    A_ineq = [A_ineq, sparse(size(A_ineq, 1), n_x)];
end
A_eq = [speye(n_x, (N+1)*n_x + N*n_u); A_dyn];
b_eq = [x0; zeros(N*n_x, 1)];

