function comparison = compare_FPK_MC_marginals_fixed()
%COMPARE_FPK_MC_MARGINALS_FIXED Compare FPK and MC marginal PDFs.
% Run with:
%   comparison = compare_FPK_MC_marginals_fixed();

% Quantitative comparison of stationary FPK and original-system MC
% marginal probability densities p(a) and p(phi).
%
% Required files in the current MATLAB folder:
%   1) selected_Omega_stationary_FPK.mat
%      produced by fpk_selected_omega_marginals.m
%   2) MC_PDF_aphi.txt
%      produced by the Fortran MC program
%
% The script:
%   * selects the requested Omega from the FPK MAT file;
%   * reconstructs the MC cell-center grids from the joint-PDF dimensions;
%   * uses cell-centered quadrature for all normalizations;
%   * interpolates the MC marginals onto the FPK grids;
%   * overlays FPK and MC p(a) and p(phi);
%   * reports peak positions, moments, high-state occupation, L1 and RMSE;
%   * estimates (but does not automatically apply) a circular phase shift.
%
% IMPORTANT:
% A nonzero best phase shift is only a diagnostic. Do not shift the MC phase
% in a publication comparison unless the two phase definitions are shown
% analytically to differ by that constant.

clc;
close all;

%% =====================================================================
% 1. USER SETTINGS
% ======================================================================

cfg.FPK_mat_file = 'selected_Omega_stationary_FPK.mat';
cfg.MC_joint_file = 'MC_PDF_aphi.txt';

cfg.Omega_target = 0.1669%0.1543;
cfg.Omega_tolerance = 5e-5;

% These limits must be identical to the Fortran histogram domain.
cfg.MC_a_min = 2.7;
cfg.MC_a_max = 3;

% Amplitude separating the low- and high-amplitude states.
% The deterministic saddle is near a=2.8--3.0 in the present case.
cfg.a_split = 2.80;

% Plot the raw theoretical phase comparison. The best circular shift is
% reported separately and is not applied by default.
cfg.apply_phase_shift = false;

% Minimum distance between the two reported peaks, in physical coordinates.
cfg.min_peak_separation_a = 0.50;
cfg.min_peak_separation_phi = 0.30;

cfg.figure_file = 'FPK_MC_marginal_comparison.png';
cfg.difference_figure_file = 'FPK_MC_density_differences.png';
cfg.save_mat_file = 'FPK_MC_comparison_results.mat';

%% =====================================================================
% 2. LOAD FPK RESULTS
% ======================================================================

if ~isfile(cfg.FPK_mat_file)
    error('Cannot find FPK file: %s',cfg.FPK_mat_file);
end

rawF = load(cfg.FPK_mat_file);

if isfield(rawF,'results')
    F = rawF.results;
else
    names = fieldnames(rawF);
    if numel(names)==1 && isstruct(rawF.(names{1}))
        F = rawF.(names{1});
    else
        error(['The MAT file does not contain a recognizable results ', ...
               'structure.']);
    end
end

required_FPK_fields = {'Omega','a','phi','pa','pphi'};
for q = 1:numel(required_FPK_fields)
    if ~isfield(F,required_FPK_fields{q})
        error('FPK results are missing field "%s".',required_FPK_fields{q});
    end
end

Omega_F_all = F.Omega(:).';
[Omega_error,kF] = min(abs(Omega_F_all-cfg.Omega_target));

if Omega_error>cfg.Omega_tolerance
    error(['Requested Omega=%.8f was not found in the FPK file. ', ...
           'Nearest available Omega is %.8f.'], ...
           cfg.Omega_target,Omega_F_all(kF));
end

Omega_used = Omega_F_all(kF);

aF = F.a(:);
phiF = F.phi(:);

paF = F.pa(:,kF);
pphiF = F.pphi(:,kF);

paF = paF(:);
pphiF = pphiF(:);

if numel(aF)~=numel(paF)
    error('FPK a and p(a) dimensions are inconsistent.');
end
if numel(phiF)~=numel(pphiF)
    error('FPK phi and p(phi) dimensions are inconsistent.');
end

daF = infer_uniform_spacing(aF,'FPK amplitude');
dphiF = infer_uniform_spacing(phiF,'FPK phase');

paF = normalize_center_density(paF,daF,'FPK p(a)');
pphiF = normalize_center_density(pphiF,dphiF,'FPK p(phi)');

%% =====================================================================
% 3. LOAD MC JOINT PDF AND COMPUTE MARGINALS
% ======================================================================

if ~isfile(cfg.MC_joint_file)
    error('Cannot find MC file: %s',cfg.MC_joint_file);
end

pM = readmatrix(cfg.MC_joint_file);

% Remove text-formatting-only empty rows or columns.
pM = pM(~all(isnan(pM),2),:);
pM = pM(:,~all(isnan(pM),1));

if isempty(pM) || any(isnan(pM),'all')
    error('MC_PDF_aphi.txt is empty or contains NaN entries.');
end

[NaM,NphiM] = size(pM);

daM = (cfg.MC_a_max-cfg.MC_a_min)/NaM;
dphiM = 2*pi/NphiM;

aM = cfg.MC_a_min + ((1:NaM)'-0.5)*daM;
phiM = ((1:NphiM)'-0.5)*dphiM;

negative_mass_MC = -sum(min(pM,0),'all')*daM*dphiM;
if negative_mass_MC>1e-12
    warning('MC joint PDF contains negative mass %.3e.',negative_mass_MC);
end
pM(pM<0) = 0;

mass_MC = sum(pM,'all')*daM*dphiM;
if ~isfinite(mass_MC) || mass_MC<=0
    error('MC joint PDF has invalid total mass.');
end
pM = pM/mass_MC;

paM = sum(pM,2)*dphiM;
pphiM = sum(pM,1)'*daM;

paM = normalize_center_density(paM,daM,'MC p(a)');
pphiM = normalize_center_density(pphiM,dphiM,'MC p(phi)');

%% =====================================================================
% 4. INTERPOLATE MC DENSITIES ONTO FPK GRIDS
% ======================================================================

% Amplitude is nonperiodic.
if numel(aM)==numel(aF) && ...
        max(abs(aM-aF))<1e-12

    paM_on_F = paM;

else

    paM_on_F = interp1( ...
        aM,paM,aF,'linear',0);

end

paM_on_F(paM_on_F<0) = 0;

paM_on_F = normalize_center_density( ...
    paM_on_F,daF,'MC p(a)');

% Phase is periodic. Extend by one period on both sides before interpolation.
if numel(phiM)==numel(phiF) && ...
        max(abs(phiM-phiF))<1e-12

    pphiM_on_F_raw = pphiM;

else

    pphiM_on_F_raw = periodic_interp( ...
        phiM,pphiM,phiF);

end

pphiM_on_F_raw(pphiM_on_F_raw<0) = 0;

pphiM_on_F_raw = normalize_center_density( ...
    pphiM_on_F_raw,dphiF,'MC p(phi)');

% Find the circular phase shift giving the smallest L1 discrepancy.
[best_shift_bins,best_shift_rad,best_shift_L1,pphiM_best] = ...
    best_circular_alignment(pphiF,pphiM_on_F_raw,dphiF);

if cfg.apply_phase_shift
    pphiM_on_F = pphiM_best;
    phase_shift_applied = best_shift_rad;
else
    pphiM_on_F = pphiM_on_F_raw;
    phase_shift_applied = 0;
end

%% =====================================================================
% 5. STATISTICS AND QUANTITATIVE ERRORS
% ======================================================================

statsF = marginal_statistics(aF,paF,phiF,pphiF, ...
    daF,dphiF,cfg.a_split, ...
    cfg.min_peak_separation_a,cfg.min_peak_separation_phi);

statsM = marginal_statistics(aF,paM_on_F,phiF,pphiM_on_F, ...
    daF,dphiF,cfg.a_split, ...
    cfg.min_peak_separation_a,cfg.min_peak_separation_phi);

L1_a = sum(abs(paF-paM_on_F))*daF;
RMSE_a = sqrt(mean((paF-paM_on_F).^2));

L1_phi = sum(abs(pphiF-pphiM_on_F))*dphiF;
RMSE_phi = sqrt(mean((pphiF-pphiM_on_F).^2));

fprintf('\n============================================================\n');
fprintf('FPK--MC MARGINAL COMPARISON\n');
fprintf('============================================================\n');
fprintf('Omega used                    = %.10f\n',Omega_used);
fprintf('FPK grid                      = %d x %d\n', ...
    numel(aF),numel(phiF));
fprintf('MC grid                       = %d x %d\n',NaM,NphiM);
fprintf('MC joint-PDF mass             = %.12f\n', ...
    sum(pM,'all')*daM*dphiM);
fprintf('Applied phase shift           = %.10f rad\n', ...
    phase_shift_applied);
fprintf('Best diagnostic phase shift   = %.10f rad (%d bins)\n', ...
    best_shift_rad,best_shift_bins);
fprintf('Best shifted phase L1         = %.6e\n',best_shift_L1);

fprintf('\nAmplitude density:\n');
fprintf('  L1 distance                 = %.6e\n',L1_a);
fprintf('  RMSE                        = %.6e\n',RMSE_a);
fprintf('  FPK mean / variance         = %.8f / %.8f\n', ...
    statsF.mean_a,statsF.var_a);
fprintf('  MC  mean / variance         = %.8f / %.8f\n', ...
    statsM.mean_a,statsM.var_a);
fprintf('  FPK P_high(a>%.3f)          = %.8f\n', ...
    cfg.a_split,statsF.P_high);
fprintf('  MC  P_high(a>%.3f)          = %.8f\n', ...
    cfg.a_split,statsM.P_high);
fprintf('  FPK amplitude peak locations=');
fprintf(' %.8f',statsF.peak_a_locations);
fprintf('\n');
fprintf('  MC  amplitude peak locations=');
fprintf(' %.8f',statsM.peak_a_locations);
fprintf('\n');

fprintf('\nPhase density:\n');
fprintf('  L1 distance                 = %.6e\n',L1_phi);
fprintf('  RMSE                        = %.6e\n',RMSE_phi);
fprintf('  FPK circular mean / R       = %.8f / %.8f\n', ...
    statsF.mean_phi_circ,statsF.R_phi);
fprintf('  MC  circular mean / R       = %.8f / %.8f\n', ...
    statsM.mean_phi_circ,statsM.R_phi);
fprintf('  FPK phase peak locations    =');
fprintf(' %.8f',statsF.peak_phi_locations);
fprintf('\n');
fprintf('  MC  phase peak locations    =');
fprintf(' %.8f',statsM.peak_phi_locations);
fprintf('\n');

%% =====================================================================
% 6. OVERLAY PLOTS
% ======================================================================

fig1 = figure('Units','centimeters','Position',[15 5 20 15]);

%tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(aF,paF,'b-','LineWidth',3.0 );
hold on;
plot(aF,paM_on_F,'ro','LineWidth',2.5);
% xline(cfg.a_split,':','LineWidth',1.0, ...
%     'DisplayName','$a_{\mathrm{s}}$');

xlabel('$a$','Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');
ylabel('$p(a)$','Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');
% title(sprintf('$p(a)$, $\\Omega=%.4f$',Omega_used), ...
%     'Interpreter','latex');
legend( {'FPK solution','MC simulation'},'FontSize',24, 'Location','northeast','Interpreter','latex','FontName','Times New Roman','Box','off');
xlim([min(aF),max(aF)]);
ylim([0,1.3*max(paF)]);
set(gca,'FontName','Times New Roman', ...
    'FontSize',24, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.5, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');
grid on;
box on;

fig2 = figure('Units','centimeters','Position',[15 5 20 15]);
plot(phiF,pphiF,'b-','LineWidth',3.0);
hold on;
plot(phiF,pphiM_on_F,'ro','LineWidth',2.5);

xlabel('$\varphi$','Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');
ylabel('$p(\varphi)$','Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');

% if cfg.apply_phase_shift
%     title(sprintf(['$p(\\varphi)$, applied MC shift $=%.3f$ rad'], ...
%         phase_shift_applied),'Interpreter','latex');
% else
%     title('$p(\varphi)$ without empirical phase shifting', ...
%         'Interpreter','latex');
% end

legend( {'FPK solution','MC simulation'}, 'FontSize',24,'Location','northeast','Interpreter','latex','FontName','Times New Roman','Box','off');
xlim([2,2.5]);
ylim([0,1.3*max(pphiF)]);
%xticks(0:pi/2:2*pi);
%xticklabels({'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'});
set(gca,'FontName','Times New Roman', ...
    'FontSize',24, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.5, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');
grid on;
box on;

% try
%     exportgraphics(fig1,cfg.figure_file,'Resolution',600);
% catch
%     print(fig1,cfg.figure_file,'-dpng','-r600');
% end

%% =====================================================================
% 7. DIFFERENCE AND CUMULATIVE-PROBABILITY PLOTS
% ======================================================================

cdfF_a = cumsum(paF)*daF;
cdfM_a = cumsum(paM_on_F)*daF;

% fig2 = figure( ...
%     'Color','w', ...
%     'Units','centimeters', ...
%     'Position',[6 6 25 10]);
% 
% tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% 
% nexttile;
% plot(aF,paM_on_F-paF,'LineWidth',3);
% yline(0,':','LineWidth',1.0);
% xlabel('$a$','Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');
% ylabel('$p_{\mathrm{MC}}(a)-p_{\mathrm{FPK}}(a)$', ...
%     'Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');
% title(sprintf('$L^1_a=%.3g$',L1_a),'Interpreter','latex');
% grid on;
% box on;
% 
% nexttile;
% plot(aF,cdfF_a,'-','LineWidth',3,'DisplayName','FPK');
% hold on;
% plot(aF,cdfM_a,'ro','LineWidth',3,'DisplayName','MC');
% xline(cfg.a_split,':','LineWidth',1.0);
% xlabel('$a$','Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');
% ylabel('Cumulative probability','Interpreter','latex', 'FontSize',24, 'FontName','Times New Roman');
% title('Amplitude cumulative distributions');
% legend('Location','best','Box','off');
% grid on;
% box on;
% 
% try
%     exportgraphics(fig2,cfg.difference_figure_file,'Resolution',600);
% catch
%     print(fig2,cfg.difference_figure_file,'-dpng','-r600');
% end

%% =====================================================================
% 8. SAVE COMPARISON DATA
% ======================================================================
% 
% comparison = struct();
% comparison.settings = cfg;
% comparison.Omega = Omega_used;
% 
% comparison.a = aF;
% comparison.phi = phiF;
% 
% comparison.pa_FPK = paF;
% comparison.pa_MC = paM_on_F;
% 
% comparison.pphi_FPK = pphiF;
% comparison.pphi_MC_raw = pphiM_on_F_raw;
% comparison.pphi_MC_compared = pphiM_on_F;
% 
% comparison.best_phase_shift_bins = best_shift_bins;
% comparison.best_phase_shift_rad = best_shift_rad;
% comparison.best_shifted_phase_L1 = best_shift_L1;
% comparison.phase_shift_applied = phase_shift_applied;
% 
% comparison.stats_FPK = statsF;
% comparison.stats_MC = statsM;
% 
% comparison.L1_a = L1_a;
% comparison.RMSE_a = RMSE_a;
% comparison.L1_phi = L1_phi;
% comparison.RMSE_phi = RMSE_phi;
% 
% save(cfg.save_mat_file,'comparison');
% 
% fprintf('\nFiles saved:\n');
% fprintf('  %s\n',cfg.figure_file);
% fprintf('  %s\n',cfg.difference_figure_file);
% fprintf('  %s\n',cfg.save_mat_file);

end

%% =====================================================================
% LOCAL FUNCTIONS
% ======================================================================

function dx = infer_uniform_spacing(x,name)
x = x(:);
if numel(x)<2
    error('%s grid must contain at least two points.',name);
end
d = diff(x);
dx = mean(d);
if max(abs(d-dx))>1e-8*max(1,abs(dx))
    error('%s grid is not uniform.',name);
end
end

function p = normalize_center_density(p,dx,name)
p = p(:);
p(p<0 & p>-1e-13) = 0;
if any(p<0)
    warning('%s contains negative values.',name);
    p(p<0) = 0;
end
mass = sum(p)*dx;
if ~isfinite(mass) || mass<=0
    error('%s has invalid mass.',name);
end
p = p/mass;
end

function pq = periodic_interp(phi,p,phi_q)
phi = phi(:);
p = p(:);
phi_q = phi_q(:);

phi_ext = [phi-2*pi;phi;phi+2*pi];
p_ext = [p;p;p];

pq = interp1(phi_ext,p_ext,phi_q,'pchip');
end

function [best_bins,best_rad,best_L1,p_best] = ...
    best_circular_alignment(p_ref,p_test,dphi)

p_ref = p_ref(:);
p_test = p_test(:);

if numel(p_ref)~=numel(p_test)
    error('Circular-alignment vectors must have equal lengths.');
end

N = numel(p_ref);
errors = zeros(N,1);

for s = 0:N-1
    p_shift = circshift(p_test,s);
    errors(s+1) = sum(abs(p_ref-p_shift))*dphi;
end

[best_L1,idx] = min(errors);
best_bins = idx-1;

% Convert shifts larger than half a period to a negative signed shift.
if best_bins>N/2
    signed_bins = best_bins-N;
else
    signed_bins = best_bins;
end

best_rad = signed_bins*dphi;
p_best = circshift(p_test,best_bins);
end

function stats = marginal_statistics(a,pa,phi,pphi, ...
    da,dphi,a_split,min_sep_a,min_sep_phi)

a = a(:);
pa = pa(:);
phi = phi(:);
pphi = pphi(:);

stats = struct();

stats.mean_a = sum(a.*pa)*da;
stats.mean_a2 = sum((a.^2).*pa)*da;
stats.var_a = max(stats.mean_a2-stats.mean_a^2,0);

stats.P_low = sum(pa(a<a_split))*da;
stats.P_high = sum(pa(a>=a_split))*da;
pnorm = stats.P_low+stats.P_high;
if pnorm>0
    stats.P_low = stats.P_low/pnorm;
    stats.P_high = stats.P_high/pnorm;
end

[idx_a,heights_a] = two_largest_separated_peaks( ...
    a,pa,min_sep_a,false);

stats.peak_a_locations = a(idx_a);
stats.peak_a_heights = heights_a;

C = sum(pphi.*cos(phi))*dphi;
S = sum(pphi.*sin(phi))*dphi;
stats.mean_phi_circ = mod(atan2(S,C),2*pi);
stats.R_phi = hypot(C,S);
stats.circular_variance_phi = 1-stats.R_phi;

[idx_phi,heights_phi] = two_largest_separated_peaks( ...
    phi,pphi,min_sep_phi,true);

stats.peak_phi_locations = phi(idx_phi);
stats.peak_phi_heights = heights_phi;
end

function [idx_selected,h_selected] = ...
    two_largest_separated_peaks(x,y,min_separation,is_periodic)

x = x(:);
y = y(:);
N = numel(y);

is_peak = false(N,1);

if is_periodic
    for i = 1:N
        im = mod(i-2,N)+1;
        ip = mod(i,N)+1;
        is_peak(i) = y(i)>y(im) && y(i)>=y(ip);
    end
else
    if N>=3
        is_peak(2:N-1) = ...
            y(2:N-1)>y(1:N-2) & y(2:N-1)>=y(3:N);
    end
end

idx_candidates = find(is_peak);

if isempty(idx_candidates)
    [~,idx] = max(y);
    idx_selected = idx;
    h_selected = y(idx);
    return;
end

[~,order] = sort(y(idx_candidates),'descend');
idx_candidates = idx_candidates(order);

idx_selected = idx_candidates(1);

for q = 2:numel(idx_candidates)
    i = idx_candidates(q);

    if is_periodic
        direct_distance = abs(x(i)-x(idx_selected(1)));
        distance = min(direct_distance,2*pi-direct_distance);
    else
        distance = abs(x(i)-x(idx_selected(1)));
    end

    if distance>=min_separation
        idx_selected(2,1) = i; %#ok<AGROW>
        break;
    end
end

idx_selected = sort(idx_selected);
h_selected = y(idx_selected);
end
