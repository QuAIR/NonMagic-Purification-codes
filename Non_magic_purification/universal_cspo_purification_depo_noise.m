clear;

%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com
% channel_magic: https://github.com/jamesrseddon/channel_magic/tree/master

%% Load A matrices for linear system %%%%
% The Amat files are from https://github.com/jamesrseddon/channel_magic/tree/master
load('Amat3.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat = A;
clear A;
[A_rows A_columns] = size(A_mat);

%% set parameter
n = 2; % the number of copies
d = 2; % dimension of qudits

pauli_array = enumeratePaulis(n+1);

purification_prob = 0.5; 
delta = 0.5; 

% depo
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

%% initial fidelity
rho = RandomDensityMatrix(d, 0, 1);
depo_rho = ApplyMap(rho, choi_depo);
init_fid = trace(depo_rho*rho)
%% cspo_primal
cvx_begin sdp quiet
    variable JE(d^(2*n),d^(2*n)) hermitian;
    % Define the optimization variable for p (as a vector)
    variable p_vec(A_columns)
    variable b_vec(A_rows)

    JN = PartialTrace(JE,3,[d^n, d, d^(n-1)]);
    f = real(trace(PartialTranspose(JN, 1, [d^n,d])*Q));

    maximize f
    subject to                                            
        JN>=0;
        PartialTrace(JN,2,[d^n,d]) <= eye(d^n);
        trace(PartialTranspose(JN, 1, [d^n,d])*R) == purification_prob;

        % Linear equality constraint with vector p
        A_mat * p_vec == calculateExpectationVec_mine(JN/d^2,pauli_array,b_vec);

        % Non-negativity of p
        p_vec >= 0;
cvx_end   


primal_fid_cspo = f/purification_prob


%% cspo_dual

cvx_begin sdp quiet
    variable Y_dual(d^(n),d^(n)) hermitian;
    variable y_vec(A_rows);
    variable x_dual;

    ob = -x_dual-trace(Y_dual)/purification_prob;
    CCC = calculateExpec(pauli_array, y_vec);

    minimize ob
    subject to     
        PartialTranspose(Q,1,[d^n,d]) + PartialTranspose(R,1, [d^n,d])*x_dual + kron(Y_dual,eye(d))- ...
         CCC <= 0;

        Y_dual == 0;

        transpose(A_mat) * y_vec <= 0;

cvx_end

dual_fid_cspo = ob


%% functions

function [dim] = get_D(n,d)
    dim = factorial(n+d-1)/(factorial(d-1)*factorial(n));
end


function [output] = calculateExpectationVec_mine(state,observables,b_vector)
[obs_rows, obs_columns, num_observables] = size(observables);

for oo = 1:num_observables
    current_obs = observables(:,:,oo);
    b_vector(oo) = trace(current_obs*state);
end

output = b_vector;

end


function [output] = calculateExpec(observables, y_vector)
    [obs_rows, obs_columns, num_observables] = size(observables);
    mat = zeros(obs_rows, obs_columns);

    for oo = 1:num_observables
        current_obs = observables(:, :, oo);
        mat = mat + y_vector(oo) * current_obs;
    end

    output = mat;
end




