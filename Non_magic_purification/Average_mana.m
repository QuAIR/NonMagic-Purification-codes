clear;
close all;

d = 3; % dimension
n = 2; % the number of copies
[psp_op1, HW1] = Generate_A(d, 1);
%% set up T state
zeta = exp(2*pi*i/9);
vt = [zeta 
      1 
      zeta^(-1)];
vt = vt / norm(vt);
rho_T = vt * vt';

%% set up H state
wh = exp(2*pi*i/3);
H = 1/sqrt(3) * [1 1 1; 1 wh wh^2; 1 wh^2 wh];
[S D] = eig(H);
v1 = S(:,1);
v2 = S(:,2);
v3 = S(:,3);
rho_H = v1*v1';

%% set up S state
v1 = [0 1 -1]';
v1 = v1 / norm(v1);
rho_S = v1 * v1';

%% set up N state
v1 = [-1 2 -1]';
v1 = v1 / norm(v1);
rho_N = v1 * v1';

n_state = 4;
rhoM(:,:,1) = rho_T;
rhoM(:,:,2) = rho_H;
rhoM(:,:,3) = rho_S;
rhoM(:,:,4) = rho_N;

delta_list = 0:0.01:1;
for epoch = 1:length(delta_list)
delta = delta_list(epoch);
depo_choi = DepolarizingChannel(d, 1 - delta);

%% Calculate average magic of the set
ave_mana = 0;
for aa = 1: length(rhoM)
ave_mana = ave_mana + Mana(psp_op1, ApplyMap(rhoM(:, :, aa), depo_choi));
end
ave_mana = ave_mana/length(rhoM); 

ave_mana_list(epoch) = ave_mana;
end

plot(delta_list, ave_mana_list)
title("Average mana")
xlabel('delta')
ylabel('Average mana')

%%
function mana = Mana(psp_op, rho)
[m, n] = size(rho);
mana = 0;
for i = 1: length(psp_op)
    wigner = trace(psp_op{i} * rho);
    mana = mana + abs(wigner)/m;
end
mana = log2(mana);
end


function [An, Tn] = Generate_A(dim, num_copy)

d = dim; % dimension
dp = 1; % used for the field Z_d for generating Heisenberg-Weyl operators
X  = GenPauli(1,0,d);
Z  = GenPauli(0,1,d);
tau = exp((d+1) * pi * 1i / d);

T1 = {};
% 1-copy Heisenberg-Weyl operators T
for a1 = 1:d
    for a2 = 1:d
        T1{end+1} = tau^(-(a1-dp) * (a2-dp)) * Z^(a1-dp) * X^(a2-dp);
    end
end

% n-copy phase-space point operators
c=1;
Tn = T1;
while c<num_copy
    Ttemp = {};
    for i=1:length(Tn)
        for j=1:d^2
            Ttemp{end+1} = Tensor(Tn{i}, T1{j});
        end
    end
    Tn = Ttemp;
c = c+1;
end

% n-copy A
An = {};
A0 = zeros(d^num_copy);
for i=1:length(Tn)
    A0 = A0 + Tn{i};
end
A0 = A0/d^num_copy;

for j=1:length(Tn)
    An{end+1} = Tn{j} * A0 * Tn{j}';
end
end