clear;
close all
%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com
% channel_magic: https://github.com/jamesrseddon/channel_magic/tree/master

%%
d = 2; % dimension
n = 2; % the number of copies


%% Load A matrices for linear system %%%%
load('Amat3.mat'); % Load A matrix for all 2-qubit stabiliser states.
Amat3 = A;
clear A;
[A_rows A_columns] = size(Amat3);
pauli_array3 = enumeratePaulis(n+1);

%%
n_state = 2;


rhoM(:, :, 1) = [1,0;0,0];
rhoM(:, :, 2) = [0.5,0.5;0.5,0.5];


p_list = 0:0.01:1;
% p_list = 0.5;
lenp = length(p_list);
purfied_f_list = [];
init_f_list = [];

for epoch = 1:lenp

clear rhoout

p = p_list(epoch);
depo_choi = DepolarizingChannel(d, 1 - p);

rhoM_depo = zeros(d, d, n_state);
rhoset = zeros(d^2, d^2, n_state);

for j=1:n_state
    rhoM_depo(:,:,j) = ApplyMap(rhoM(:,:,j), depo_choi);
    rhoset(:,:,j) = kron(rhoM_depo(:,:,j), rhoM_depo(:,:,j));
end

initialf = 0;
for j=1:n_state
    initialf = initialf + trace(rhoM_depo(:,:,j) * rhoM(:,:,j));
end

p_distill = 0.6;


cvx_begin sdp quiet
    variable JN(d^3, d^3) hermitian
    % Define the optimization variable for p (as a vector)
    variable p_vec(A_columns)
    variable b_vec(A_rows)

    f = 0;
    init_prob = 0;

    for j=1:n_state
        rhoout(:,:,j) = PartialTrace(JN * kron(rhoset(:,:,j).', eye(d)), 1, [d^2 d]);
        f = f + real(trace(rhoout(:,:,j)/p_distill * rhoM(:,:,j)));
    end

    maximise f
        JN >= 0;
        PartialTrace(JN, 2, [d^2 d]) <= eye(d^2);

        for j=1:n_state
            init_prob = init_prob + trace(rhoout(:,:,j));
        end
        init_prob/n_state == p_distill;

        % Linear equality constraint with vector p
        Amat3 * p_vec == calculateExpectationVec_mine(JN/d^2,pauli_array3,b_vec);

        % Non-negativity of p
        p_vec>= 0;

cvx_end

purfied_f = f/n_state
init_f = initialf/n_state
purfied_f_list(epoch) = purfied_f;
init_f_list(epoch) = init_f;


end

figure;
plot(p_list, purfied_f_list,'r', 'DisplayName',sprintf('CSPO purifying', p_distill));
hold on;
plot(p_list, init_f_list, 'b', 'DisplayName','no purifying');

title('Purifying fidelity vs. error rate', 'Interpreter', 'latex');
xlabel('$\delta$', 'Interpreter', 'latex');
ylabel(sprintf('Fidelity(successful probability p = %.1f )', p_distill), 'Interpreter', 'latex');
legend show;
grid on;


%% functions


function [output] = calculateExpectationVec_mine(state,observables,b_vector)
[obs_rows, obs_columns, num_observables] = size(observables);

for oo = 1:num_observables
    current_obs = observables(:,:,oo);
    b_vector(oo) = trace(current_obs*state);
end

output = b_vector;

end


