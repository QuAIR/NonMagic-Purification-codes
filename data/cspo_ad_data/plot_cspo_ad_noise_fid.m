clear;

load fid_ad_noise_p_05_1.mat

gamma = 0:0.005:1;

init_f_list = 1/3 + (1 + sqrt(1 - gamma)).^2/6;

color_p05 = [0.4660, 0.6740, 0.1880]; 
color_p10 = [0.8500, 0.3250, 0.0980]; 
color_base = [0.4, 0.4, 0.4]; 

% figure p = 0.5
figure;

ax = gca; 

lw = 2; % LineWidth
mk_idx = 1:5:length(gamma); 
ms = 7; % MarkerSize

plot(gamma, f_cspo_05, '-', 'Color', color_p05, 'LineWidth', lw, ...
    'Marker', 'o', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', '$\mathcal{A}_{\mathrm{CSPO}} (p=0.5)$');
hold on

plot(gamma, f_cpwp_05, '--', 'Color', color_p05, 'LineWidth', lw, ...
    'Marker', 'x', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', 'CPTN $(p=0.5)$');
hold on

plot(gamma, init_f_list, ':', 'Color', color_base, 'LineWidth', 3, ...
    'DisplayName', 'No purify');
hold on

xlabel('$\gamma$', 'Interpreter', 'latex', 'Interpreter', 'latex');
ylabel('Fidelity', 'Interpreter', 'latex', 'Interpreter', 'latex');
title('$p=0.5$', 'Interpreter', 'latex');

lgd = legend('Interpreter', 'latex', 'Location', 'southwest');
set(gca, 'FontName', 'Times');

% ========= FONT SIZE CONTROL =========
ax.FontSize = 20;         % ticks on both x and y
% Legend text
lgd.FontSize = 20;

hold on;
box on;
grid on;

print(gcf, '-dsvg', '-vector', 'fid_cspo_ad_05.svg');


% figure p = 1
figure; 
plot(gamma, f_cspo_1, '-', 'Color', color_p10, 'LineWidth', lw, ...
    'Marker', 'o', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', '$\mathcal{A}_{\mathrm{CSPO}} (p=1.0)$');
hold on
plot(gamma, f_cpwp_1, '--', 'Color', color_p10, 'LineWidth', lw, ...
    'Marker', 'x', 'MarkerIndices', mk_idx, 'MarkerSize', ms, ...
    'DisplayName', 'CPTN $(p=1.0)$');
hold on

plot(gamma, init_f_list, ':', 'Color', color_base, 'LineWidth', 3, ...
    'DisplayName', 'No purify');
hold on

xlabel('$\gamma$', 'Interpreter', 'latex', 'Interpreter', 'latex');
ylabel('Fidelity', 'Interpreter', 'latex', 'Interpreter', 'latex');
title('$p=1$', 'Interpreter', 'latex');

lgd = legend('Interpreter', 'latex', 'Location', 'southwest');
set(gca, 'FontName', 'Times');

% ========= FONT SIZE CONTROL =========
ax.FontSize = 20;         % ticks on both x and y
% Legend text
lgd.FontSize = 20;

hold on;
box on;
grid on;

print(gcf, '-dsvg', '-vector', 'fid_cspo_ad_1.svg');



