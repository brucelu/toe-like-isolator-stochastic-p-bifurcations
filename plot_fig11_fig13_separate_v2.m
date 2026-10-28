function results = plot_fig11_fig13_separate_v2(output_dir)
%PLOT_FIG11_FIG13_DETERMINISTIC_MATCHED
% Self-contained code for the revised Fig. 11 and Fig. 13.
%
% The code avoids using global stationary branch weights in the weak-noise
% bistable regime. Instead, it shows:
%
% Fig. 11
%   (a) deterministic low-amplitude, saddle, and high-amplitude slow-flow
%       fixed-point branches versus xi2 under matched effective linear
%       damping;
%   (b) the amplitude-dependent effective damping rate
%
%           Gamma(a;xi2) = rho0 + xi2*b2*a^2/4.
%
% Fig. 13
%   (a) deterministic mono-/bistability map in the (Omega,xi1) plane with
%       xi2 fixed;
%   (b) maximum stable-response amplitude in the same plane;
%   (c) deterministic mono-/bistability map in the (Omega,xi2) plane under
%       matched effective linear damping;
%   (d) maximum stable-response amplitude in the matched-damping plane.
%
% Usage
% -----
% results = plot_fig11_fig13_deterministic_matched;
%
% or
%
% results = plot_fig11_fig13_deterministic_matched( ...
%     'Fig11_Fig13_deterministic_results');
%
% All figures are saved as JPG, PNG, PDF, and editable FIG files.
%
% This program requires only standard MATLAB functions.

clc;
close all;

if nargin < 1 || isempty(output_dir)
    output_dir = 'Fig11_Fig13_separate_results_v2';
end

validateattributes(output_dir,{'char','string'},{'nonempty'}, ...
    mfilename,'output_dir',1);
output_dir = char(output_dir);

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% 1. SLOW-FLOW PARAMETERS
par.k1 = 0.0313;
par.k2 = -0.0609;
par.k3 = 0.0381;
par.b0 = 0.5000;
par.b2 = 0.5417;
par.epsv = 0.1;
par.Zhat = 0.001;

par.wn = sqrt(par.k1);
par.Zs = par.Zhat/(par.epsv^3);
par.kappa = ...
    (3*par.k3-10*par.k2^2/(3*par.wn^2))/(2*par.wn);

% Reference damping pair used to define the matched condition.
par.xi1_ref = 0.005;
par.xi2_ref = 0.005;
par.rhohat_ref = par.xi1_ref + par.xi2_ref*par.b0;
par.rhos_target = par.rhohat_ref/(par.epsv^2);

fprintf('\n====================================================\n');
fprintf('REVISED FIG. 11 AND SEPARATE FIG. 13 PANELS\n');
fprintf('====================================================\n');
fprintf('wn                  = %.12f\n',par.wn);
fprintf('kappa               = %.12f\n',par.kappa);
fprintf('reference xi1       = %.6f\n',par.xi1_ref);
fprintf('reference xi2       = %.6f\n',par.xi2_ref);
fprintf('matched rho_s       = %.12f\n',par.rhos_target);
fprintf('matched relation    : xi1 = %.6f - %.6f*xi2\n', ...
    par.epsv^2*par.rhos_target,par.b0);

%% 2. PLOT SETTINGS
cfg.line_width = 2.6;
cfg.marker_size = 6;
cfg.font_size = 24;
cfg.label_font_size = 24;
cfg.figure_width_cm = 32;
cfg.figure_height_cm = 13.5;

% Deterministic response-state classification:
% 1 = low-amplitude monostability
% 2 = bistability
% 3 = high-amplitude monostability
cfg.branch_threshold = 2.5;
cfg.state_colors = [ ...
    0.20,0.45,0.85; ...
    0.88,0.25,0.20; ...
    0.95,0.70,0.18];
cfg.state_ticks = [1,2,3];
cfg.state_tick_labels = { ...
    'Low-amplitude monostable', ...
    'Bistable', ...
    'High-amplitude monostable'};

%% 3. FIG. 11: MATCHED-DAMPING NONLINEAR-DAMPING MECHANISM
Omega_fig11 = 0.1530;
xi2_fig11 = linspace(0.001,0.010,181);

xi1_fig11 = ...
    par.epsv^2*par.rhos_target-par.b0*xi2_fig11;

if any(xi1_fig11<0)
    error('The selected xi2 range produces a negative matched xi1.');
end

% Determine an amplitude threshold from the reference fixed points. This
% threshold is used only if a parameter value has a single stable root.
par_ref = par;
par_ref.xi1 = par.xi1_ref;
par_ref.xi2 = par.xi2_ref;
par_ref.rhos = par.rhos_target;
ref_info = deterministic_fixed_points(Omega_fig11,par_ref);
ref_stable = sort(ref_info.a_roots(ref_info.stable));

if numel(ref_stable)>=2
    branch_threshold = 0.5*(ref_stable(1)+ref_stable(end));
else
    branch_threshold = 2.5;
end

a_low = NaN(size(xi2_fig11));
a_saddle = NaN(size(xi2_fig11));
a_high = NaN(size(xi2_fig11));
n_stable_fig11 = zeros(size(xi2_fig11));

for k = 1:numel(xi2_fig11)
    local_par = par;
    local_par.xi2 = xi2_fig11(k);
    local_par.xi1 = xi1_fig11(k);
    local_par.rhos = ...
        (local_par.xi1+local_par.xi2*local_par.b0)/(local_par.epsv^2);

    info = deterministic_fixed_points(Omega_fig11,local_par);
    [a_low(k),a_saddle(k),a_high(k),n_stable_fig11(k)] = ...
        split_branches(info,branch_threshold);
end

% Relative nonlinear-damping contribution evaluated on each stable branch:
%
% R_nl^(j) = xi2*b2*(a_j^*)^2/(4*rho_s),  j=L,H.
%
% The result is reported as a percentage.
Rnl_low = 100*xi2_fig11.*par.b2.*a_low.^2 ...
    /(4*par.rhos_target);
Rnl_high = 100*xi2_fig11.*par.b2.*a_high.^2 ...
    /(4*par.rhos_target);

% ----------------------------------------------------------------------
% Fig. 11(a): deterministic branches
% ----------------------------------------------------------------------
fig11a = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 22 20]);

hold on;

h_low = plot(xi2_fig11,a_low, ...
    'LineWidth',cfg.line_width);
h_high = plot(xi2_fig11,a_high, ...
    'LineWidth',cfg.line_width);
h_saddle = plot(xi2_fig11,a_saddle,'--', ...
    'LineWidth',cfg.line_width);

xlabel('$\xi_2$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
ylabel('Fixed-point amplitude $a^*$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);

legend([h_low,h_high,h_saddle], ...
    {'Low-amplitude stable branch', ...
     'High-amplitude stable branch', ...
     'Saddle branch'}, ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',24, ...
    'Box','off');
% 坐标轴统一设置
ax = gca;
set(ax, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'TickDir', 'in');

xlim([min(xi2_fig11),max(xi2_fig11)]);
grid on;
box on;
%format_axes(gca,cfg);

save_all_formats(fig11a,output_dir, ...
    'Fig11a_deterministic_branches');

% ----------------------------------------------------------------------
% Fig. 11(b): branch-resolved nonlinear-damping contribution
% ----------------------------------------------------------------------
fig11b = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 22 20]);

hold on;

h_Rlow = plot(xi2_fig11,Rnl_low, ...
    'LineWidth',cfg.line_width);
h_Rhigh = plot(xi2_fig11,Rnl_high, ...
    'LineWidth',cfg.line_width);

xlabel('$\xi_2$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
ylabel('Relative nonlinear-damping contribution (\%)', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);

legend([h_Rlow,h_Rhigh], ...
    {'Low-amplitude stable branch', ...
     'High-amplitude stable branch'}, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'FontSize',24, ...
    'Box','off');
% 坐标轴统一设置
ax = gca;
set(ax, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'TickDir', 'in');

xlim([min(xi2_fig11),max(xi2_fig11)]);
grid on;
box on;
%format_axes(gca,cfg);

save_all_formats(fig11b,output_dir, ...
    'Fig11b_branch_resolved_nonlinear_damping');

% Save Fig. 11 data.
T11 = table( ...
    xi2_fig11(:), ...
    xi1_fig11(:), ...
    par.rhos_target*ones(numel(xi2_fig11),1), ...
    a_low(:), ...
    a_saddle(:), ...
    a_high(:), ...
    n_stable_fig11(:), ...
    Rnl_low(:), ...
    Rnl_high(:), ...
    'VariableNames', ...
    {'xi2','matched_xi1','rho_s','a_low','a_saddle','a_high', ...
     'n_stable','Rnl_low_percent','Rnl_high_percent'});

writetable(T11, ...
    fullfile(output_dir,'Fig11_matched_branch_data.csv'));

%% 4. FIG. 13: FOUR SEPARATE DETERMINISTIC PARAMETER-PLANE FIGURES
Omega_grid = linspace(0.135,0.165,241);
xi1_grid = linspace(0.001,0.010,181);
xi2_grid = linspace(0.001,0.010,181);

xi2_fixed = 0.005;

% (Omega,xi1) plane.
state_xi1 = NaN(numel(xi1_grid),numel(Omega_grid));
a_upper_xi1 = NaN(size(state_xi1));

fprintf('\nComputing (Omega,xi1) plane: %d x %d points...\n', ...
    numel(xi1_grid),numel(Omega_grid));

for iy = 1:numel(xi1_grid)
    local_par = par;
    local_par.xi1 = xi1_grid(iy);
    local_par.xi2 = xi2_fixed;
    local_par.rhos = ...
        (local_par.xi1+local_par.xi2*local_par.b0) ...
        /(local_par.epsv^2);

    for ix = 1:numel(Omega_grid)
        info = deterministic_fixed_points(Omega_grid(ix),local_par);

        [state_code,upper_amp] = classify_response_state( ...
            info,cfg.branch_threshold);

        state_xi1(iy,ix) = state_code;
        a_upper_xi1(iy,ix) = upper_amp;
    end
end

% Matched (Omega,xi2) plane.
matched_xi1_grid = ...
    par.epsv^2*par.rhos_target-par.b0*xi2_grid;

valid_xi2 = matched_xi1_grid>=0;

state_xi2 = NaN(numel(xi2_grid),numel(Omega_grid));
a_upper_xi2 = NaN(size(state_xi2));

fprintf('Computing matched (Omega,xi2) plane: %d x %d points...\n', ...
    numel(xi2_grid),numel(Omega_grid));

for iy = 1:numel(xi2_grid)
    if ~valid_xi2(iy)
        continue;
    end

    local_par = par;
    local_par.xi2 = xi2_grid(iy);
    local_par.xi1 = matched_xi1_grid(iy);
    local_par.rhos = par.rhos_target;

    for ix = 1:numel(Omega_grid)
        info = deterministic_fixed_points(Omega_grid(ix),local_par);

        [state_code,upper_amp] = classify_response_state( ...
            info,cfg.branch_threshold);

        state_xi2(iy,ix) = state_code;
        a_upper_xi2(iy,ix) = upper_amp;
    end
end

% ----------------------------------------------------------------------
% Fig. 13(a): response-state map in the (Omega,xi1) plane
% ----------------------------------------------------------------------
fig13a = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 28 22]);
%     'MenuBar','none', ...
%     'ToolBar','none'

ax13a = axes(fig13a, ...
    'Units','normalized', ...
    'Position',[0.10 0.14 0.69 0.80]);

imagesc(ax13a,Omega_grid,1e3*xi1_grid,state_xi1);
set(ax13a,'YDir','normal');
colormap(ax13a,cfg.state_colors);
caxis(ax13a,[0.5,3.5]);

cb = colorbar(ax13a,'eastoutside');
cb.Units = 'normalized';
cb.Position = [0.82 0.14 0.035 0.80];
cb.Ticks = cfg.state_ticks;
cb.TickLabels = {'Low','Bistability','High'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 24;
cb.Title.String = 'State';
cb.Title.Interpreter = 'latex';
cb.Title.FontSize = 24;

hold(ax13a,'on');
contour(ax13a,Omega_grid,1e3*xi1_grid,state_xi1, ...
    [1.5,1.5],'k-','LineWidth',2);
contour(ax13a,Omega_grid,1e3*xi1_grid,state_xi1, ...
    [2.5,2.5],'k-','LineWidth',2);

xlabel(ax13a,'$\Omega$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
ylabel(ax13a,'$10^3\xi_1$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
% 坐标轴统一设置
ax = gca;
set(ax, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'TickDir', 'in','YTick', [2, 4, 6, 8, 10],'xTick', [0.135, 0.145, 0.155, 0.165]);

%format_axes(ax13a,cfg);

save_all_formats(fig13a,output_dir, ...
    'Fig13a_Omega_xi1_response_states');

% ----------------------------------------------------------------------
% Fig. 13(b): upper stable-branch amplitude in the (Omega,xi1) plane
% ----------------------------------------------------------------------
fig13b = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 28 22]);
% 'MenuBar','none', ...
%     'ToolBar','none'

ax13b = axes(fig13b, ...
    'Units','normalized', ...
    'Position',[0.10 0.14 0.73 0.80]);

imagesc(ax13b,Omega_grid,1e3*xi1_grid,a_upper_xi1);
set(ax13b,'YDir','normal');
colormap(ax13b,parula(256));

cb = colorbar(ax13b,'eastoutside');
cb.Units = 'normalized';
cb.Position = [0.86 0.14 0.030 0.80];
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 24;
cb.Title.String = '$a_{\mathrm{upper}}^*$';
cb.Title.Interpreter = 'latex';
cb.Title.FontSize = 24;

hold(ax13b,'on');
contour(ax13b,Omega_grid,1e3*xi1_grid,state_xi1, ...
    [1.5,1.5],'k-','LineWidth',2);
contour(ax13b,Omega_grid,1e3*xi1_grid,state_xi1, ...
    [2.5,2.5],'k-','LineWidth',2);

xlabel(ax13b,'$\Omega$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
ylabel(ax13b,'$10^3\xi_1$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);

% 坐标轴统一设置
ax = gca;
set(ax, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'TickDir', 'in','YTick', [2, 4, 6, 8, 10],'xTick', [0.135, 0.145, 0.155, 0.165]);
%format_axes(ax13b,cfg);

save_all_formats(fig13b,output_dir, ...
    'Fig13b_Omega_xi1_upper_stable_amplitude');

% ----------------------------------------------------------------------
% Fig. 13(c): matched response-state map in the (Omega,xi2) plane
% ----------------------------------------------------------------------
fig13c = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 28 22]);
%     'MenuBar','none', ...
%     'ToolBar','none');

ax13c = axes(fig13c, ...
    'Units','normalized', ...
    'Position',[0.10 0.14 0.69 0.80]);

imagesc(ax13c,Omega_grid,1e3*xi2_grid,state_xi2);
set(ax13c,'YDir','normal');
colormap(ax13c,cfg.state_colors);
caxis(ax13c,[0.5,3.5]);

cb = colorbar(ax13c,'eastoutside');
cb.Units = 'normalized';
cb.Position = [0.82 0.14 0.035 0.80];
cb.Ticks = cfg.state_ticks;
cb.TickLabels = {'Low','Bistability','High'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 24;
cb.Title.String = 'State';
cb.Title.Interpreter = 'latex';
cb.Title.FontSize = 24;

hold(ax13c,'on');
contour(ax13c,Omega_grid,1e3*xi2_grid,state_xi2, ...
    [1.5,1.5],'k-','LineWidth',2);
contour(ax13c,Omega_grid,1e3*xi2_grid,state_xi2, ...
    [2.5,2.5],'k-','LineWidth',2);

xlabel(ax13c,'$\Omega$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
ylabel(ax13c,'$10^3\xi_2$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
% 坐标轴统一设置
ax = gca;
set(ax, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'TickDir', 'in','YTick', [2, 4, 6, 8, 10],'xTick', [0.135, 0.145, 0.155, 0.165]);

%format_axes(ax13c,cfg);

save_all_formats(fig13c,output_dir, ...
    'Fig13c_Omega_xi2_matched_response_states');

% ----------------------------------------------------------------------
% Fig. 13(d): matched upper stable-branch amplitude
% ----------------------------------------------------------------------
fig13d = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 28 22]);
%  'MenuBar','none', ...
%     'ToolBar','none'

ax13d = axes(fig13d, ...
    'Units','normalized', ...
    'Position',[0.10 0.14 0.73 0.80]);

imagesc(ax13d,Omega_grid,1e3*xi2_grid,a_upper_xi2);
set(ax13d,'YDir','normal');
colormap(ax13d,parula(256));

cb = colorbar(ax13d,'eastoutside');
cb.Units = 'normalized';
cb.Position = [0.86 0.14 0.030 0.80];
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 24;
cb.Title.String = '$a_{\mathrm{upper}}^*$';
cb.Title.Interpreter = 'latex';
cb.Title.FontSize = 24;

hold(ax13d,'on');
contour(ax13d,Omega_grid,1e3*xi2_grid,state_xi2, ...
    [1.5,1.5],'k-','LineWidth',2);
contour(ax13d,Omega_grid,1e3*xi2_grid,state_xi2, ...
    [2.5,2.5],'k-','LineWidth',2);

xlabel(ax13d,'$\Omega$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
ylabel(ax13d,'$10^3\xi_2$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.label_font_size);
% 坐标轴统一设置
ax = gca;
set(ax, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'TickDir', 'in','YTick', [2, 4, 6, 8, 10],'xTick', [0.135, 0.145, 0.155, 0.165]);

%format_axes(ax13d,cfg);

save_all_formats(fig13d,output_dir, ...
    'Fig13d_Omega_xi2_matched_upper_stable_amplitude');

%% 5. SAVE NUMERICAL DATA
results = struct();
results.parameters = par;
results.settings = cfg;

results.fig11.Omega = Omega_fig11;
results.fig11.xi2 = xi2_fig11;
results.fig11.xi1 = xi1_fig11;
results.fig11.a_low = a_low;
results.fig11.a_saddle = a_saddle;
results.fig11.a_high = a_high;
results.fig11.n_stable = n_stable_fig11;
results.fig11.Rnl_low_percent = Rnl_low;
results.fig11.Rnl_high_percent = Rnl_high;

results.fig13.Omega = Omega_grid;
results.fig13.xi1 = xi1_grid;
results.fig13.xi2 = xi2_grid;
results.fig13.matched_xi1 = matched_xi1_grid;
results.fig13.state_xi1 = state_xi1;
results.fig13.a_upper_xi1 = a_upper_xi1;
results.fig13.state_xi2 = state_xi2;
results.fig13.a_upper_xi2 = a_upper_xi2;

save(fullfile(output_dir,'Fig11_Fig13_separate_data.mat'), ...
    'results','-v7.3');

[OMEGA_XI1,XI1_MESH] = meshgrid(Omega_grid,xi1_grid);

T13_xi1 = table( ...
    OMEGA_XI1(:), ...
    XI1_MESH(:), ...
    xi2_fixed*ones(numel(OMEGA_XI1),1), ...
    state_xi1(:), ...
    a_upper_xi1(:), ...
    'VariableNames', ...
    {'Omega','xi1','xi2_fixed','state_code','a_upper_stable'});

writetable(T13_xi1, ...
    fullfile(output_dir,'Fig13_Omega_xi1_data.csv'));

[OMEGA_XI2,XI2_MESH] = meshgrid(Omega_grid,xi2_grid);
XI1_MATCH_MESH = repmat(matched_xi1_grid(:),1,numel(Omega_grid));

T13_xi2 = table( ...
    OMEGA_XI2(:), ...
    XI2_MESH(:), ...
    XI1_MATCH_MESH(:), ...
    par.rhos_target*ones(numel(OMEGA_XI2),1), ...
    state_xi2(:), ...
    a_upper_xi2(:), ...
    'VariableNames', ...
    {'Omega','xi2','matched_xi1','rho_s', ...
     'state_code','a_upper_stable'});

writetable(T13_xi2, ...
    fullfile(output_dir,'Fig13_Omega_xi2_matched_data.csv'));

fprintf('\n====================================================\n');
fprintf('CALCULATION COMPLETED\n');
fprintf('====================================================\n');
fprintf('Output folder:\n%s\n',output_dir);
fprintf('\nSeparate figures:\n');
fprintf('  Fig11a_deterministic_branches.*\n');
fprintf('  Fig11b_branch_resolved_nonlinear_damping.*\n');
fprintf('  Fig13a_Omega_xi1_response_states.*\n');
fprintf('  Fig13b_Omega_xi1_upper_stable_amplitude.*\n');
fprintf('  Fig13c_Omega_xi2_matched_response_states.*\n');
fprintf('  Fig13d_Omega_xi2_matched_upper_stable_amplitude.*\n');
fprintf('\nState codes:\n');
fprintf('  1 = low-amplitude monostability\n');
fprintf('  2 = bistability\n');
fprintf('  3 = high-amplitude monostability\n');
end

%% =====================================================================
% Deterministic slow-flow fixed points
% =====================================================================
function info = deterministic_fixed_points(Omega,par)

sigma = (Omega-par.wn)/(par.epsv^2);
c = par.xi2*par.b2/4;
d = par.kappa/4;
F = par.Zs/(2*par.wn);

% Polynomial in y=a^2:
%
% (c^2+d^2)y^3 + 2(rho*c-sigma*d)y^2
% + (rho^2+sigma^2)y - F^2 = 0.
coef = [ ...
    c^2+d^2, ...
    2*(par.rhos*c-sigma*d), ...
    par.rhos^2+sigma^2, ...
    -F^2];

yroot = roots(coef);

imag_tol = 1e-8;
positive_tol = 1e-12;
mask = abs(imag(yroot)) <= ...
    imag_tol.*max(1,abs(real(yroot))) & ...
    real(yroot)>positive_tol;

a_roots = sort(sqrt(real(yroot(mask))));
a_roots = remove_close_values(a_roots,1e-7);

stable = false(size(a_roots));
phi_roots = NaN(size(a_roots));
max_real_eig = NaN(size(a_roots));

for j = 1:numel(a_roots)
    aj = a_roots(j);

    sin_phi = aj*(par.rhos+c*aj^2)/F;
    cos_phi = aj*(-sigma+d*aj^2)/F;

    radius = hypot(sin_phi,cos_phi);
    if radius>0
        sin_phi = sin_phi/radius;
        cos_phi = cos_phi/radius;
    end

    phi_j = mod(atan2(sin_phi,cos_phi),2*pi);

    J11 = -par.rhos-3*c*aj^2;
    J12 = F*cos(phi_j);
    J21 = -2*d*aj-F*cos(phi_j)/(aj^2);
    J22 = -(F/aj)*sin(phi_j);

    eigJ = eig([J11,J12;J21,J22]);
    max_real_eig(j) = max(real(eigJ));
    stable(j) = max_real_eig(j)<-1e-9;
    phi_roots(j) = phi_j;
end

info = struct();
info.a_roots = a_roots;
info.phi_roots = phi_roots;
info.stable = stable;
info.max_real_eig = max_real_eig;
end

%% =====================================================================
% Separate low, saddle, and high branches
% =====================================================================
function [a_low,a_saddle,a_high,n_stable] = ...
    split_branches(info,branch_threshold)

a_low = NaN;
a_saddle = NaN;
a_high = NaN;

stable_amp = sort(info.a_roots(info.stable));
unstable_amp = sort(info.a_roots(~info.stable));
n_stable = numel(stable_amp);

if numel(stable_amp)>=2
    a_low = stable_amp(1);
    a_high = stable_amp(end);

    between = unstable_amp( ...
        unstable_amp>a_low & unstable_amp<a_high);

    if ~isempty(between)
        midpoint = 0.5*(a_low+a_high);
        [~,idx] = min(abs(between-midpoint));
        a_saddle = between(idx);
    elseif ~isempty(unstable_amp)
        a_saddle = unstable_amp(1);
    end
elseif numel(stable_amp)==1
    if stable_amp(1)<branch_threshold
        a_low = stable_amp(1);
    else
        a_high = stable_amp(1);
    end

    if ~isempty(unstable_amp)
        a_saddle = unstable_amp(1);
    end
end
end

%% =====================================================================
% Classify deterministic response state
% =====================================================================
function [state_code,upper_amp] = ...
    classify_response_state(info,branch_threshold)

stable_amp = sort(info.a_roots(info.stable));

state_code = NaN;
upper_amp = NaN;

if isempty(stable_amp)
    return;
end

upper_amp = stable_amp(end);

if numel(stable_amp)>=2
    state_code = 2;
elseif stable_amp(1)<branch_threshold
    state_code = 1;
else
    state_code = 3;
end
end

%% =====================================================================
% Remove nearly repeated roots
% =====================================================================
function x_unique = remove_close_values(x,tol)

x = sort(x(:).');

if isempty(x)
    x_unique = x;
    return;
end

x_unique = x(1);

for k = 2:numel(x)
    scale = max([1,abs(x(k)),abs(x_unique(end))]);
    if abs(x(k)-x_unique(end))>tol*scale
        x_unique(end+1) = x(k); %#ok<AGROW>
    end
end
end

%% =====================================================================
% Axes formatting
% =====================================================================
% function format_axes(ax,cfg)
% 
% 
% ax = gca;
% set(ax, ...
%     'FontSize', 24, ...
%     'FontName', 'Times New Roman', ...
%     'LineWidth', 1.5, ...
%     'Box', 'on', ...
%     'Layer', 'top', ...
%     'TickDir', 'in');
% set(ax, ...
%     'FontName','Times New Roman', ...
%     'FontSize',cfg.font_size, ...
%     'TickLabelInterpreter','latex', ...
%     'LineWidth',1.5, ...
%     'TickDir','in', ...
%     'Layer','top', ...
%     'Box','on');
% end

%% =====================================================================
% Save an individual axes panel
% =====================================================================
function save_axes_formats(ax_handle,output_dir,base_name)

jpg_file = fullfile(output_dir,[base_name,'.jpg']);
png_file = fullfile(output_dir,[base_name,'.png']);
pdf_file = fullfile(output_dir,[base_name,'.pdf']);

try
    exportgraphics(ax_handle,jpg_file,'Resolution',600);
catch
    warning('Could not export individual JPG panel: %s',base_name);
end

try
    exportgraphics(ax_handle,png_file,'Resolution',600);
catch
    warning('Could not export individual PNG panel: %s',base_name);
end

try
    exportgraphics(ax_handle,pdf_file,'ContentType','vector');
catch
    warning('Could not export individual PDF panel: %s',base_name);
end
end

%% =====================================================================
% Save figures
% =====================================================================
function save_all_formats(fig_handle,output_dir,base_name)

jpg_file = fullfile(output_dir,[base_name,'.jpg']);
png_file = fullfile(output_dir,[base_name,'.png']);
pdf_file = fullfile(output_dir,[base_name,'.pdf']);
fig_file = fullfile(output_dir,[base_name,'.fig']);

try
    exportgraphics(fig_handle,jpg_file,'Resolution',600);
catch
    print(fig_handle,jpg_file,'-djpeg','-r600');
end

try
    exportgraphics(fig_handle,png_file,'Resolution',600);
catch
    print(fig_handle,png_file,'-dpng','-r600');
end

try
    exportgraphics(fig_handle,pdf_file,'ContentType','vector');
catch
    print(fig_handle,pdf_file,'-dpdf','-painters');
end

savefig(fig_handle,fig_file);
end
