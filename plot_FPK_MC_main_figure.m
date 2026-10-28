function metrics = plot_FPK_MC_main_figure()
%PLOT_FPK_MC_MAIN_FIGURE
% Publication-ready comparison between the stationary FPK solution and the
% original-system Monte Carlo simulation.
%
% Main figure:
%   (a) marginal amplitude density p(a);
%   (b) marginal phase density p(phi);
%   (c) 50%, 80%, and 95% highest-density-region contours in the (a,phi)
%       plane.
%
% Required files:
%   selected_Omega_stationary_FPK.mat
%   MC_PDF_aphi.txt
%
% The FPK MAT file must contain results.Omega, results.a, results.phi,
% and a joint-PDF field, preferably results.P.
%
% Run:
%   metrics = plot_FPK_MC_main_figure();

clc;
close all;

%% 1. USER SETTINGS
cfg.FPK_mat_file = 'selected_Omega_stationary_FPK.mat';
cfg.MC_joint_file = 'MC_PDF_aphi.txt';

cfg.Omega_target = 0.1669;
cfg.Omega_tolerance = 5.0e-5;

% Must match the Fortran histogram domain.
cfg.MC_a_min = 2.7;
cfg.MC_a_max = 3.0;

% Plot windows.
cfg.a_plot_limits = [2.7, 3.0];
cfg.phi_plot_limits = [2.0, 2.5];

% Keep zero if phase correction was already applied in Fortran.
% If raw MC phase is loaded and phi_corrected = phi_raw-delta_phi,
% set this value to delta_phi.
cfg.additional_MC_phase_shift = 0.0;

% Highest-density regions enclosing these probability masses.
cfg.HDR_probability_masses = [0.50, 0.80, 0.95];

% Plot every nth MC point. Set [] for automatic selection.
cfg.marker_stride = [];

% Publication formatting.
cfg.figure_width_cm = 18.5;
cfg.figure_height_cm = 6.8;
cfg.font_name = 'Times New Roman';
cfg.font_size = 8.5;
cfg.axis_line_width = 0.8;
cfg.curve_line_width = 1.6;
cfg.marker_size = 3.6;

% Output files.
cfg.output_base = 'FPK_MC_main_validation_figure';
cfg.metrics_file = 'FPK_MC_main_validation_metrics.mat';

%% 2. LOAD FPK SOLUTION
if ~isfile(cfg.FPK_mat_file)
    error('Cannot find the FPK file: %s', cfg.FPK_mat_file);
end

rawF = load(cfg.FPK_mat_file);
if isfield(rawF, 'results')
    F = rawF.results;
else
    names = fieldnames(rawF);
    if numel(names) == 1 && isstruct(rawF.(names{1}))
        F = rawF.(names{1});
    else
        error('The MAT file does not contain a recognizable results structure.');
    end
end

required_fields = {'Omega', 'a', 'phi'};
for q = 1:numel(required_fields)
    if ~isfield(F, required_fields{q})
        error('The FPK results are missing field "%s".', required_fields{q});
    end
end

Omega_all = F.Omega(:).';
[Omega_error, kOmega] = min(abs(Omega_all - cfg.Omega_target));
if Omega_error > cfg.Omega_tolerance
    error(['Requested Omega=%.8f was not found. ', ...
           'The nearest available value is %.8f.'], ...
           cfg.Omega_target, Omega_all(kOmega));
end
Omega_used = Omega_all(kOmega);

aF = F.a(:);
phiF = F.phi(:);
NaF = numel(aF);
NphiF = numel(phiF);

daF = uniform_spacing(aF, 'FPK amplitude');
dphiF = uniform_spacing(phiF, 'FPK phase');

pF = extract_fpk_joint_pdf(F, kOmega, NaF, NphiF);
pF = normalize_joint_pdf(pF, daF, dphiF, 'FPK joint PDF');

% Derive all marginals from the same joint PDF used in panel (c).
paF = sum(pF, 2) * dphiF;
pphiF = sum(pF, 1).' * daF;
paF = normalize_density_1d(paF, daF, 'FPK p(a)');
pphiF = normalize_density_1d(pphiF, dphiF, 'FPK p(phi)');

%% 3. LOAD MC JOINT DENSITY
if ~isfile(cfg.MC_joint_file)
    error('Cannot find the MC file: %s', cfg.MC_joint_file);
end

pM = readmatrix(cfg.MC_joint_file);
pM = pM(~all(isnan(pM), 2), :);
pM = pM(:, ~all(isnan(pM), 1));

if isempty(pM) || any(isnan(pM), 'all')
    error('MC_PDF_aphi.txt is empty or contains NaN entries.');
end

[NaM, NphiM] = size(pM);
daM = (cfg.MC_a_max - cfg.MC_a_min) / NaM;
dphiM = 2*pi / NphiM;

aM = cfg.MC_a_min + ((1:NaM).' - 0.5) * daM;
phiM = ((1:NphiM).' - 0.5) * dphiM;
pM = normalize_joint_pdf(pM, daM, dphiM, 'MC joint PDF');

%% 4. OPTIONAL PHASE-COORDINATE CORRECTION
if abs(cfg.additional_MC_phase_shift) > 0
    pM = shift_joint_pdf_in_phase( ...
        aM, phiM, pM, cfg.additional_MC_phase_shift);
    pM = normalize_joint_pdf( ...
        pM, daM, dphiM, 'phase-corrected MC joint PDF');
end

%% 5. PLACE MC DENSITY ON FPK GRID
same_grid = ...
    NaM == NaF && ...
    NphiM == NphiF && ...
    max(abs(aM - aF)) < 1.0e-12 && ...
    max(abs(phiM - phiF)) < 1.0e-12;

if same_grid
    pM_on_F = pM;
else
    pM_on_F = interpolate_joint_pdf_periodic( ...
        aM, phiM, pM, aF, phiF);
    pM_on_F(pM_on_F < 0) = 0;
    pM_on_F = normalize_joint_pdf( ...
        pM_on_F, daF, dphiF, 'interpolated MC joint PDF');
end

paM = sum(pM_on_F, 2) * dphiF;
pphiM = sum(pM_on_F, 1).' * daF;
paM = normalize_density_1d(paM, daF, 'MC p(a)');
pphiM = normalize_density_1d(pphiM, dphiF, 'MC p(phi)');

%% 6. QUANTITATIVE COMPARISON
L1_a = sum(abs(paM - paF)) * daF;
L1_phi = sum(abs(pphiM - pphiF)) * dphiF;
RMSE_a = sqrt(mean((paM - paF).^2));
RMSE_phi = sqrt(mean((pphiM - pphiF).^2));

L1_joint = sum(abs(pM_on_F - pF), 'all') * daF * dphiF;
RMSE_joint = sqrt(mean((pM_on_F - pF).^2, 'all'));
hellinger_joint = sqrt( ...
    0.5 * sum((sqrt(pM_on_F) - sqrt(pF)).^2, 'all') ...
    * daF * dphiF);

inner_product = sum(pM_on_F .* pF, 'all') * daF * dphiF;
norm_F = sqrt(sum(pF.^2, 'all') * daF * dphiF);
norm_M = sqrt(sum(pM_on_F.^2, 'all') * daF * dphiF);
density_correlation = inner_product / (norm_F * norm_M);

mean_a_F = sum(aF .* paF) * daF;
mean_a_M = sum(aF .* paM) * daF;
var_a_F = sum((aF - mean_a_F).^2 .* paF) * daF;
var_a_M = sum((aF - mean_a_M).^2 .* paM) * daF;

[~, ~, mean_phi_F, R_phi_F] = circular_statistics(phiF, pphiF, dphiF);
[~, ~, mean_phi_M, R_phi_M] = circular_statistics(phiF, pphiM, dphiF);

[~, idx_mode_F] = max(pF(:));
[iF, jF] = ind2sub(size(pF), idx_mode_F);
[~, idx_mode_M] = max(pM_on_F(:));
[iM, jM] = ind2sub(size(pM_on_F), idx_mode_M);
mode_F = [aF(iF), phiF(jF)];
mode_M = [aF(iM), phiF(jM)];

masses = cfg.HDR_probability_masses(:).';
levels_F = highest_density_thresholds(pF, daF, dphiF, masses);
levels_M = highest_density_thresholds(pM_on_F, daF, dphiF, masses);

fprintf('\n============================================================\n');
fprintf('PUBLICATION FIGURE: FPK--MC VALIDATION\n');
fprintf('============================================================\n');
fprintf('Omega                         = %.10f\n', Omega_used);
fprintf('FPK grid                      = %d x %d\n', NaF, NphiF);
fprintf('MC grid                       = %d x %d\n', NaM, NphiM);

fprintf('\nMarginal amplitude density:\n');
fprintf('  L1 distance                 = %.6e\n', L1_a);
fprintf('  RMSE                        = %.6e\n', RMSE_a);
fprintf('  FPK mean / variance         = %.8f / %.8e\n', mean_a_F, var_a_F);
fprintf('  MC  mean / variance         = %.8f / %.8e\n', mean_a_M, var_a_M);

fprintf('\nMarginal phase density:\n');
fprintf('  L1 distance                 = %.6e\n', L1_phi);
fprintf('  RMSE                        = %.6e\n', RMSE_phi);
fprintf('  FPK circular mean / R       = %.8f / %.8f\n', mean_phi_F, R_phi_F);
fprintf('  MC  circular mean / R       = %.8f / %.8f\n', mean_phi_M, R_phi_M);

fprintf('\nJoint density:\n');
fprintf('  L1 distance                 = %.6e\n', L1_joint);
fprintf('  RMSE                        = %.6e\n', RMSE_joint);
fprintf('  Hellinger distance          = %.6e\n', hellinger_joint);
fprintf('  Density correlation         = %.8f\n', density_correlation);
fprintf('  FPK mode (a,phi)            = (%.8f, %.8f)\n', mode_F(1), mode_F(2));
fprintf('  MC  mode (a,phi)            = (%.8f, %.8f)\n', mode_M(1), mode_M(2));

%% 7. PUBLICATION-READY THREE-PANEL FIGURE
color_FPK = [0.0000, 0.3000, 0.7500];
color_MC = [0.8500, 0.2000, 0.1000];

if isempty(cfg.marker_stride)
    cfg.marker_stride = max(1, round(numel(aF) / 45));
end
idx_a_marker = 1:cfg.marker_stride:numel(aF);
idx_phi_marker = 1:cfg.marker_stride:numel(phiF);

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2.0, 2.0, cfg.figure_width_cm, cfg.figure_height_cm], ...
    'Renderer', 'painters');

tl = tiledlayout(fig, 1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Panel (a): p(a)
ax1 = nexttile(tl, 1);
plot(ax1, aF, paF, '-', ...
    'Color', color_FPK, ...
    'LineWidth', cfg.curve_line_width, ...
    'DisplayName', 'FPK solution');
hold(ax1, 'on');
plot(ax1, aF(idx_a_marker), paM(idx_a_marker), 'o', ...
    'Color', color_MC, ...
    'MarkerEdgeColor', color_MC, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', cfg.marker_size, ...
    'LineWidth', 0.9, ...
    'DisplayName', 'MC simulation');
xlabel(ax1, '$a$', 'Interpreter', 'latex');
ylabel(ax1, '$p(a)$', 'Interpreter', 'latex');
title(ax1, '(a) Amplitude density', ...
    'Interpreter', 'latex', 'FontWeight', 'normal');
xlim(ax1, cfg.a_plot_limits);
ylim(ax1, [0, 1.08 * max([paF; paM])]);
legend(ax1, 'Location', 'northwest', ...
    'Interpreter', 'latex', 'Box', 'off', ...
    'FontSize', cfg.font_size - 0.5);
grid(ax1, 'on');
box(ax1, 'on');

% Panel (b): p(phi)
ax2 = nexttile(tl, 2);
plot(ax2, phiF, pphiF, '-', ...
    'Color', color_FPK, ...
    'LineWidth', cfg.curve_line_width);
hold(ax2, 'on');
plot(ax2, phiF(idx_phi_marker), pphiM(idx_phi_marker), 'o', ...
    'Color', color_MC, ...
    'MarkerEdgeColor', color_MC, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', cfg.marker_size, ...
    'LineWidth', 0.9);
xlabel(ax2, '$\phi$', 'Interpreter', 'latex');
ylabel(ax2, '$p(\phi)$', 'Interpreter', 'latex');
title(ax2, '(b) Phase density', ...
    'Interpreter', 'latex', 'FontWeight', 'normal');
xlim(ax2, cfg.phi_plot_limits);
ylim(ax2, [0, 1.08 * max([pphiF; pphiM])]);
grid(ax2, 'on');
box(ax2, 'on');

% Panel (c): highest-density-region contours
ax3 = nexttile(tl, 3);
hold(ax3, 'on');

% masses = [0.50,0.80,0.95]. The outer 95% contour is drawn first.
draw_order = numel(masses):-1:1;
contour_line_widths = [2.0, 1.5, 1.0];

for q = draw_order
    contour(ax3, phiF, aF, pF, ...
        [levels_F(q), levels_F(q)], ...
        'LineStyle', '-', ...
        'LineWidth', contour_line_widths(q), ...
        'LineColor', color_FPK, ...
        'HandleVisibility', 'off');

    contour(ax3, phiF, aF, pM_on_F, ...
        [levels_M(q), levels_M(q)], ...
        'LineStyle', '--', ...
        'LineWidth', contour_line_widths(q), ...
        'LineColor', color_MC, ...
        'HandleVisibility', 'off');
end

% Clean legend: no data1/data2 entries.
plot(ax3, nan, nan, '-', ...
    'Color', color_FPK, ...
    'LineWidth', cfg.curve_line_width, ...
    'DisplayName', 'FPK solution');
plot(ax3, nan, nan, '--', ...
    'Color', color_MC, ...
    'LineWidth', cfg.curve_line_width, ...
    'DisplayName', 'MC simulation');

xlabel(ax3, '$\phi$', 'Interpreter', 'latex');
ylabel(ax3, '$a$', 'Interpreter', 'latex');
title(ax3, '(c) Joint-density contours', ...
    'Interpreter', 'latex', 'FontWeight', 'normal');
xlim(ax3, cfg.phi_plot_limits);
ylim(ax3, cfg.a_plot_limits);
legend(ax3, 'Location', 'northeast', ...
    'Interpreter', 'latex', 'Box', 'off', ...
    'FontSize', cfg.font_size - 0.5);
text(ax3, 0.03, 0.04, ...
    'Outer to inner: 95\%, 80\%, 50\%', ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', cfg.font_name, ...
    'FontSize', cfg.font_size - 1.0, ...
    'VerticalAlignment', 'bottom');
grid(ax3, 'on');
box(ax3, 'on');

set([ax1, ax2, ax3], ...
    'FontName', cfg.font_name, ...
    'FontSize', cfg.font_size, ...
    'LineWidth', cfg.axis_line_width, ...
    'TickLabelInterpreter', 'latex', ...
    'Layer', 'top');

title(tl, sprintf('$\\Omega=%.4f$', Omega_used), ...
    'Interpreter', 'latex', ...
    'FontName', cfg.font_name, ...
    'FontSize', cfg.font_size + 1.0, ...
    'FontWeight', 'normal');

%% 8. EXPORT
pdf_file = [cfg.output_base, '.pdf'];
png_file = [cfg.output_base, '.png'];
fig_file = [cfg.output_base, '.fig'];

try
    exportgraphics(fig, pdf_file, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');
catch
    print(fig, pdf_file, '-dpdf', '-painters');
end

try
    exportgraphics(fig, png_file, ...
        'Resolution', 600, ...
        'BackgroundColor', 'white');
catch
    print(fig, png_file, '-dpng', '-r600');
end

savefig(fig, fig_file);

%% 9. SAVE METRICS
metrics = struct();
metrics.settings = cfg;
metrics.Omega = Omega_used;
metrics.a = aF;
metrics.phi = phiF;
metrics.p_FPK = pF;
metrics.p_MC = pM_on_F;
metrics.pa_FPK = paF;
metrics.pa_MC = paM;
metrics.pphi_FPK = pphiF;
metrics.pphi_MC = pphiM;
metrics.L1_a = L1_a;
metrics.RMSE_a = RMSE_a;
metrics.L1_phi = L1_phi;
metrics.RMSE_phi = RMSE_phi;
metrics.L1_joint = L1_joint;
metrics.RMSE_joint = RMSE_joint;
metrics.Hellinger_joint = hellinger_joint;
metrics.density_correlation = density_correlation;
metrics.mean_a_FPK = mean_a_F;
metrics.mean_a_MC = mean_a_M;
metrics.var_a_FPK = var_a_F;
metrics.var_a_MC = var_a_M;
metrics.mean_phi_FPK = mean_phi_F;
metrics.mean_phi_MC = mean_phi_M;
metrics.R_phi_FPK = R_phi_F;
metrics.R_phi_MC = R_phi_M;
metrics.mode_FPK = mode_F;
metrics.mode_MC = mode_M;
metrics.HDR_probability_masses = masses;
metrics.HDR_levels_FPK = levels_F;
metrics.HDR_levels_MC = levels_M;
save(cfg.metrics_file, 'metrics');

fprintf('\nFiles saved:\n');
fprintf('  %s\n', pdf_file);
fprintf('  %s\n', png_file);
fprintf('  %s\n', fig_file);
fprintf('  %s\n', cfg.metrics_file);

end

%% LOCAL FUNCTIONS
function p = extract_fpk_joint_pdf(F, kOmega, Na, Nphi)
candidate_fields = {'P','p','joint_pdf','pdf_joint','Pjoint','stationary_pdf'};
field_found = '';
for q = 1:numel(candidate_fields)
    if isfield(F, candidate_fields{q})
        field_found = candidate_fields{q};
        break;
    end
end
if isempty(field_found)
    error(['The FPK MAT file contains no stationary joint PDF. ', ...
        'Before saving the FPK results, add "results.P = P_all;".']);
end
raw = F.(field_found);
if iscell(raw)
    if kOmega > numel(raw)
        error('The FPK joint-PDF cell array does not contain this Omega.');
    end
    p = raw{kOmega};
else
    p = select_numeric_joint_slice(raw, kOmega, Na, Nphi);
end
if isequal(size(p), [Na, Nphi])
elseif isequal(size(p), [Nphi, Na])
    p = p.';
else
    error('Extracted FPK joint PDF has size %dx%d; expected %dx%d.', ...
        size(p,1), size(p,2), Na, Nphi);
end
end

function p = select_numeric_joint_slice(raw, kOmega, Na, Nphi)
dims = size(raw);
if ismatrix(raw)
    p = raw;
    return;
end
if ndims(raw) ~= 3
    error('Unsupported FPK joint-PDF array with %d dimensions.', ndims(raw));
end
if dims(1) == Na && dims(2) == Nphi
    p = raw(:,:,kOmega);
elseif dims(1) == Nphi && dims(2) == Na
    p = raw(:,:,kOmega).';
elseif dims(2) == Na && dims(3) == Nphi
    p = squeeze(raw(kOmega,:,:));
elseif dims(2) == Nphi && dims(3) == Na
    p = squeeze(raw(kOmega,:,:)).';
else
    error('Cannot identify the FPK joint-PDF dimensions [%s].', num2str(dims));
end
end

function p = normalize_joint_pdf(p, da, dphi, name)
p = double(p);
p(p < 0 & p > -1.0e-13) = 0;
negative_mass = -sum(min(p,0), 'all') * da * dphi;
if negative_mass > 1.0e-10
    warning('%s contains negative mass %.3e.', name, negative_mass);
end
p(p < 0) = 0;
mass = sum(p, 'all') * da * dphi;
if ~isfinite(mass) || mass <= 0
    error('%s has invalid total probability mass.', name);
end
p = p / mass;
end

function p = normalize_density_1d(p, dx, name)
p = double(p(:));
p(p < 0 & p > -1.0e-13) = 0;
if any(p < 0)
    warning('%s contains negative values; they were clipped to zero.', name);
    p(p < 0) = 0;
end
mass = sum(p) * dx;
if ~isfinite(mass) || mass <= 0
    error('%s has invalid total probability mass.', name);
end
p = p / mass;
end

function dx = uniform_spacing(x, name)
x = x(:);
if numel(x) < 2
    error('%s grid must contain at least two points.', name);
end
d = diff(x);
dx = mean(d);
if max(abs(d - dx)) > 1.0e-8 * max(1, abs(dx))
    error('%s grid is not uniform.', name);
end
end

function pq = interpolate_joint_pdf_periodic(a, phi, p, a_q, phi_q)
a = a(:);
phi = phi(:);
a_q = a_q(:);
phi_q = phi_q(:);
phi_ext = [phi - 2*pi; phi; phi + 2*pi];
p_ext = [p, p, p];
[PHI_ext, A_ext] = meshgrid(phi_ext, a);
[PHI_q, A_q] = meshgrid(phi_q, a_q);
pq = interp2(PHI_ext, A_ext, p_ext, PHI_q, A_q, 'linear', 0);
end

function p_shifted = shift_joint_pdf_in_phase(a, phi, p, delta_phi)
% If phi_corrected = phi_raw-delta_phi, then
% p_corrected(phi) = p_raw(phi+delta_phi).
a = a(:);
phi = phi(:);
phi_query = mod(phi + delta_phi, 2*pi);
phi_ext = [phi - 2*pi; phi; phi + 2*pi];
p_ext = [p, p, p];
[PHI_ext, A_ext] = meshgrid(phi_ext, a);
[PHI_q, A_q] = meshgrid(phi_query, a);
p_shifted = interp2(PHI_ext, A_ext, p_ext, PHI_q, A_q, 'linear', 0);
end

function levels = highest_density_thresholds(P, da, dphi, masses)
masses = masses(:).';
if any(masses <= 0) || any(masses >= 1)
    error('Highest-density-region masses must lie strictly between 0 and 1.');
end
values = sort(P(:), 'descend');
cumulative_mass = cumsum(values) * da * dphi;
levels = zeros(size(masses));
for q = 1:numel(masses)
    idx = find(cumulative_mass >= masses(q), 1, 'first');
    if isempty(idx)
        levels(q) = values(end);
    else
        levels(q) = values(idx);
    end
end
end

function [C, S, mean_phi, R] = circular_statistics(phi, pphi, dphi)
phi = phi(:);
pphi = pphi(:);
C = sum(pphi .* cos(phi)) * dphi;
S = sum(pphi .* sin(phi)) * dphi;
mean_phi = mod(atan2(S, C), 2*pi);
R = hypot(C, S);
end
