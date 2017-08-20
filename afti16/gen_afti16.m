% AFTI-16 2x2 system

% (C) 2003 by Alberto Bemporad
% edited Jan 30, 2017 by Lorenzo Stella

n_x = 4;
n_u = 2;

% Continuous-time model

Ac = [-.0151 -60.5651 0 -32.174;
      -.0001 -1.3411 .9929 0;
       .00018 43.2541 -.86939 0;
       0      0       1      0];
Bc = [-2.516 -13.136;
     -.1689 -.2514;
     -17.251 -1.5766;
     0        0];
Cc = [0 1 0 0;
     0 0 0 1];
Dc = [0 0;
     0 0];

sys = ss(Ac,Bc,Cc,Dc);

% Discrete-time model

Ts = .05; % Sampling time

model = c2d(sys, Ts);

A = model.A;
B = model.B;
C = model.C;
D = model.D;
