clear;

load fid_qutrit_p_01_06.mat

load ave_mana_depo_qutrit.mat
ave_mana = max(ave_mana_list, 0);

p_list = 0:0.01:1;
init_f_list = 1 - 2/3*p_list;   % initial fidelity

figure;

% Colors for each axis/curve
leftColor  = [0, 0.4470, 0.7410];   % blue-ish
rightColor = [0.8500, 0.3250, 0.0980];   % reddish

% --- Left y-axis: purification gap ---
yyaxis left
ax = gca;                   % current axes handle
ax.YColor = leftColor;      % set left y-axis color

p1 = plot(p_list, f_cpwp_01 - init_f_list, ':', ...
    'LineWidth', 2.5, ...
    'Color', leftColor, ...
    'DisplayName', 'Fidelity gap');
ylabel('Fidelity gap $F^{\mathcal{A}}_{\mathcal{D}_\delta} - \lambda_0$', 'Interpreter', 'latex');

ylim([-0.01 0.14]);

% --- Right y-axis: average mana ---
yyaxis right
ax.YColor = rightColor;     % set right y-axis color

p2 = plot(p_list, ave_mana, '-', ...
    'LineWidth', 2.5, ...
    'Color', rightColor, ...
    'DisplayName', 'Average mana');
ylabel('Average mana', 'Interpreter', 'latex');

% optional: customize right y-axis limits
ylim([-0.05 0.75]);

% --- Common settings ---
xlabel('$\delta$', 'Interpreter', 'latex');
title('$p=0.1$', ...
      'Interpreter', 'latex');

% Legend (you can still customize its location)
lgd = legend([p1 p2], 'Interpreter', 'latex', 'Location', 'none');
lgd.Position = [0.38 0.76 0.18 0.08];  % adjust as you like

% ========= FONT SIZE CONTROL =========
ax.FontSize = 20;         % ticks on both x and y
% Legend text
lgd.FontSize = 20;

set(gca, 'FontName', 'Times');
grid on;

print(gcf, '-dsvg', '-vector', 'fid_gap_ave_mana_cpwp.svg');