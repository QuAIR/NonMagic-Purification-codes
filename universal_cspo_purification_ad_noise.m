clear;

%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com
% channel_magic: https://github.com/jamesrseddon/channel_magic/tree/master

%% Load A matrices for linear system %%%%
load('Amat3.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat = A;
clear A;
[A_rows A_columns] = size(A_mat);

%% set parameter
n = 2; % the number of copies
d = 2; % dimension of qudits

pauli_array = enumeratePaulis(n+1);

purification_prob = 1; 
puri_fid_record = [];

deltas = 0:0.005:1;
for i = 1:length(deltas)
delta = deltas(i); 

%% depo or ad
% choi_depo = DepolarizingChannel(d, 1 - delta);
AD = amplitudeDamping(delta);
K1 = AD(:, :, 1);
K2 = AD(:, :, 2);
E1 = ChoiMatrix({K1; K2});
choi_ad = E1;

JN = kron(choi_ad, choi_ad);
JN = PermuteSystems(JN, [1,3,2,4], [d,d,d,d]);

JI = MaxEntangled(d,0,0)*MaxEntangled(d,0,0)';
J = PermuteSystems(kron(JN,JI), [1,2,5,3,4,6], [d,d,d,d,d,d]);
PI = full(SymmetricProjection(d,n+1));
Q = ApplyMap(PI,J) / get_D(n+1,d);

PI = full(SymmetricProjection(d,n));
R0 = ApplyMap(PI, JN) / get_D(n,d);
R = kron(R0,eye(d));


%% cspo_primal
cvx_begin sdp quiet
    variable JE(d^(2*n),d^(2*n)) hermitian;
    % Define the optimization variable for p (as a vector)
    variable p_vec(A_columns)
    variable b_vec(A_rows)

    JN = PartialTrace(JE,3,[d^n, d, d^(n-1)]);% 3qubit
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

puri_fid_record(i) = primal_fid_cspo;
end
figure;
plot(deltas, puri_fid_record,'DisplayName','purify')
hold on 
plot(deltas, 1/3 + (1 + sqrt(1 - deltas)).^2/6,'DisplayName','no purify')
legend;

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



