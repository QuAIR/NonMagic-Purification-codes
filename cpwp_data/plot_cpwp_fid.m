clear;

load fid_qutrit_p_01_06.mat

p_list = 0:0.01:1;

init_f_list = 1 - 2/3*p_list;


color_p01 = [0, 0.4470, 0.7410]; 
color_p06 = [0.4660, 0.6740, 0.1880]; 
color_base = [0.4, 0.4, 0.4]; 


% figure p = 0.1
figure;
ax = gca;

lw = 2; % LineWidth
mk_idx = 1:3:length(p_list); 
ms = 7; % MarkerSize

plot(p_list, f_cpwp_01, '-', 'Color', color_p01, 'LineWidth', lw, ...
    'Marker', 'o', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', '$\mathcal{A}_{\mathrm{CPWP}} (p=0.1)$');
hold on

plot(p_list, f_cptn_01, '--', 'Color', color_p01, 'LineWidth', lw, ...
    'Marker', 'x', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', 'CPTN $(p=0.1)$');
hold on
plot(p_list, init_f_list, ':', 'LineWidth', 2.5,'Color', color_base, 'DisplayName','No purify');

title('$p=0.1$', 'Interpreter', 'latex');
xlabel('$\delta$', 'Interpreter', 'latex');
ylabel('Fidelity', 'Interpreter', 'latex');
lgd = legend('Location', 'southwest', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times');
set(gca, 'YLim', [0.3 1.05]);

% ========= FONT SIZE CONTROL =========
ax.FontSize = 20;         % ticks on both x and y
% Legend text
lgd.FontSize = 20;

grid on;
print(gcf, '-dsvg', '-vector', 'fid_cpwp_01.svg');



% figure p = 0.6
figure;
ax = gca;

lw = 2; % LineWidth
mk_idx = 1:3:length(p_list); 
ms = 7; % MarkerSize
plot(p_list, f_cpwp_06, '-', 'Color', color_p06, 'LineWidth', lw, ...
    'Marker', 'o', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', '$\mathcal{A}_{\mathrm{CPWP}} (p=0.6)$');
hold on
plot(p_list, f_cptn_06, '--', 'Color', color_p06, 'LineWidth', lw, ...
    'Marker', 'x', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', 'CPTN $(p=0.6)$');
hold on


plot(p_list, init_f_list, ':', 'LineWidth', 2.5,'Color', color_base, 'DisplayName','No purify');


title('$p=0.6$', 'Interpreter', 'latex');
xlabel('$\delta$', 'Interpreter', 'latex');
ylabel('Fidelity', 'Interpreter', 'latex');
lgd = legend('Location', 'southwest', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times');
set(gca, 'YLim', [0.3 1.05]);

% ========= FONT SIZE CONTROL =========
ax.FontSize = 20;         % ticks on both x and y
% Legend text
lgd.FontSize = 20;

grid on;

print(gcf, '-dsvg', '-vector', 'fid_cpwp_06.svg');



