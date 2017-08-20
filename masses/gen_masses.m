function [A, B, ts] = gen_masses(n_controls)

if mod(n_controls,2) ~= 0
    error('number of actuators must be even');
end

k = 1;          % spring constant 
lam = 0;        % damping constant 
a = -2*k;
b = -2*lam;
c = k;
d = lam;

n_masses = n_controls*2;

n = 2*n_masses; % state dimension
m = n_controls; % input dimension

Acts = [zeros(n/2), eye(n/2); ...
    spdiags([c*ones(n/2,1), a*ones(n/2,1), c*ones(n/2,1)], [-1, 0, 1], n/2, n/2), ...
    spdiags([d*ones(n/2,1), b*ones(n/2,1), d*ones(n/2,1)], [-1, 0, 1], n/2, n/2)
];

Bcts = [zeros(n_masses,m); ...
    kron(speye(n_controls/2), [1 0; 0 1; -1 0; 0 -1]);
];

% convert to discrete-time system
ts = 0.5;       % sampling time
A = expm(ts*Acts);
B = (Acts\(expm(ts*Acts)-eye(n)))*Bcts;
