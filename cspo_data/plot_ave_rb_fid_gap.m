clear;


load fid_qubit_p_01_06.mat

ave_rb = zeros(1, length(f_cspo_01));

p_list = 0:0.01:1;
init_f_list = 1 - 1/2*p_list;   % initial fidelity

figure;

% Colors for each axis/curve
leftColor  = [0, 0.4470, 0.7410];   % blue-ish
rightColor = [0.8500, 0.3250, 0.0980];   % reddish

% --- Left y-axis: purification gap ---
yyaxis left
ax = gca;                   % current axes handle
ax.YColor = leftColor;      % set left y-axis color

p1 = plot(p_list, f_cspo_01 -  init_f_list, ':', ...
    'LineWidth', 2.5, ...
    'Color', leftColor, ...
    'DisplayName', 'Fidelity gap');

y1l = ylabel('Fidelity gap $F^{\mathcal{A}}_{\mathcal{D}_\delta} - \lambda_0$', 'Interpreter', 'latex');

% optional: customize left y-axis limits
ylim([-0.01 0.26]);

% --- Right y-axis: average mana ---
yyaxis right
ax.YColor = rightColor;     % set right y-axis color

p2 = plot(p_list, ave_rb, '-', ...
    'LineWidth', 2.5, ...
    'Color', rightColor, ...
    'DisplayName', 'Average robustness');
y2l = ylabel('Average log robustness', 'Interpreter', 'latex');

% optional: customize right y-axis limits
ylim([-0.05 0.75]);

% --- Common settings ---
xl = xlabel('$\delta$', 'Interpreter', 'latex');


tl = title('$p=0.1$', 'Interpreter', 'latex');

% Legend (you can still customize its location)
lgd = legend([p1 p2], 'Interpreter', 'latex', 'Location', 'none');
lgd.Position = [0.33 0.77 0.18 0.08];  % adjust as you like

% ========= FONT SIZE CONTROL =========
ax.FontSize = 20;         % ticks on both x and y
% Legend text
lgd.FontSize = 20;


set(gca, 'FontName', 'Times');
grid on;

print(gcf, '-dsvg', '-vector', 'fid_gap_ave_rb_cspo.svg');