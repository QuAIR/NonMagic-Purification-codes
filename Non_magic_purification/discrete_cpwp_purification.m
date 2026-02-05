clear;
close all;
%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com

%%
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

%%

v_0 = [1,0,0];
v_1 = [0,1,0];
v_2 = [0,0,1];
rho_0 = v_0.'*v_0;

v_plus = (v_0 + v_1 + v_2)/sqrt(3);
rho_plus = v_plus.'*v_plus;

%%
n_state = 4;
rhoM(:,:,1) = rho_T;
rhoM(:,:,2) = rho_H;
rhoM(:,:,3) = rho_S;
rhoM(:,:,4) = rho_N;

% p_list = 0:0.01:1;
p_list = 0.5;
purfied_f_list = [];
init_f_list = [];

p_distill = 0.6;

for epoch = 1:length(p_list)
%% 
clear rhoout
clear psp_coffi

depo_choi = DepolarizingChannel(d, 1 - p);


for j=1:n_state
    rhoM_depo(:,:,j) = ApplyMap(rhoM(:,:,j), depo_choi);
    rhoset(:,:,j) = kron(rhoM_depo(:,:,j), rhoM_depo(:,:,j));
end

initialf = 0;
for j=1:n_state
    initialf = initialf + trace(rhoM_depo(:,:,j) * rhoM(:,:,j));
end


cvx_begin sdp quiet
cvx_solver sedumi
    variable JN(d^3, d^3) hermitian

    f = 0;
    for j=1:n_state
        rhoout(:,:,j) = PartialTrace(JN * kron(rhoset(:,:,j).', eye(d)), 1, [d^2 d]);
        f = f + real(trace(rhoout(:,:,j)/p_distill * rhoM(:,:,j)));
    end
    maximise f
        JN >= 0;
        PartialTrace(JN, 2, [d^2 d]) <= eye(d^2);

    for j = 1:length(psp_op1)
        for k = 1:length(psp_op1)
            for i = 1: length(psp_op1)
            t = Tensor(psp_op1{j}.' , psp_op1{k}.', psp_op1{i});
            psp_coffi(j,k,i) = real(trace(t * JN));
            psp_coffi(j,k,i) >= 0;            
            end
        end
    end

        for j=1:n_state
            trace(rhoout(:,:,j)) == p_distill;
        end
cvx_end

purfied_f = f/n_state
init_f = initialf/n_state
purfied_f_list(epoch) = purfied_f;
init_f_list(epoch) = init_f;

end

figure;

plot(p_list, purfied_f_list,'r', 'DisplayName',sprintf('CPWP purifying', p_distill));
hold on;
plot(p_list, init_f_list, 'b', 'DisplayName','no purifying');




title('Purifying fidelity vs. error rate', 'Interpreter', 'latex');
xlabel('$\delta$', 'Interpreter', 'latex');
ylabel(sprintf('Fidelity(successful probability p = %.1f )', p_distill), 'Interpreter', 'latex');
legend show;
grid on;

%%
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
