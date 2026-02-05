clear;
%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com

%% set parameter
n = 2; % the number of copies
d = 3; % dimension of qudits
[psp_op1, HW1] = Generate_A(d, 1);

% there are d^(2*n) psp operators, with each dimension d^n
[psp_op_n, HWn] = Generate_A(d, n);

delta = 0.5;

choi_depo = DepolarizingChannel(d, 1 - delta);
JN = choi_depo;
for i = 1:n-1
    JN = kron(JN, choi_depo);
end

JN = PermuteSystems(JN, [2*[0:n-1] + 1, 2*[0:n-1] + 2], [d*ones(1,2*n)]);

JI = MaxEntangled(d,0,0)*MaxEntangled(d,0,0)';
J = PermuteSystems(kron(JN,JI), [1:n,2*n+1,n+1:2*n,2*n+2]);
PI = full(SymmetricProjection(d,n+1));
Q = ApplyMap(PI,J) / get_D(n+1,d);

PI = full(SymmetricProjection(d,n));
R0 = ApplyMap(PI, JN) / get_D(n,d);
R = kron(R0,eye(d));

purification_prob = 0.5;

rho = RandomDensityMatrix(d, 0, 1);
depo_rho = ApplyMap(rho, choi_depo);
init_fid = trace(depo_rho*rho)

%% primal SDP
% %%% cpwp %%%
cvx_begin sdp quiet

    variable JN(d^(n+1),d^(n+1)) hermitian;

    f = real(trace(PartialTranspose(JN, 1, [d^n,d])*Q));

    maximize f
    subject to                                            
        JN>=0;
        PartialTrace(JN,2,[d^n,d]) <=eye(d^n);
        trace(PartialTranspose(JN, 1, [d^n,d])*R) == purification_prob;

        % % % % % % CPWP map
        % there are d^(2*(n+2)) constrain, with each dimension d^(n+1)
        for j = 1:length(psp_op_n)
            for k = 1:length(psp_op1)           
                t = kron(psp_op_n{j}.', psp_op1{k});
                real(trace(t * JN)) >= 0;
            end
        end


cvx_end

primal_fid = f/purification_prob
% 
% 
% 
%% dual SDP
clear JN
clear t

t_uv = 0;
cvx_begin sdp quiet

variable Y(d^n, d^n) hermitian
variable x 
variable c(length(psp_op_n), length(psp_op1)) 

f_dual = -x - trace(Y)/purification_prob;

pt_Q = PartialTranspose(Q, 1, [d^n,d]);
pt_R = PartialTranspose(R, 1, [d^n,d]) * x;
kr_Y = kron(Y, eye(d));

for j = 1:length(psp_op_n)
    for k = 1:length(psp_op1)
        t_uv = t_uv + c(j,k)*kron(psp_op_n{j}.', psp_op1{k})/d;
    end
end

minimize  f_dual

subject to 
Y == 0;


for j = 1:length(psp_op_n)
    for k = 1:length(psp_op1)
            c(j,k) <= 0;
    end
end

pt_Q + pt_R + kr_Y - t_uv <= 0;

cvx_end

f_dual




%% functions

function [dim] = get_D(n,d)
    dim = factorial(n+d-1)/(factorial(d-1)*factorial(n));
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