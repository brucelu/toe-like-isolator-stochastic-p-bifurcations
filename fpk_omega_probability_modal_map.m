function results = fpk_omega_probability_modal_map(Omega_scan,S_eta)
%FPK_OMEGA_PROBABILITY_MODAL_MAP
% Construct the stationary amplitude-probability modal map in the
% (Omega,a) plane and identify the unimodal-bimodal P-bifurcation boundaries.
%
% The figure contains:
%   (a) column-normalized stationary amplitude density p(a;Omega), with all
%       detected modal amplitudes and P-bifurcation boundaries superposed;
%   (b) detected modal-amplitude branches. Blue and red denote unimodal and
%       bimodal stationary amplitude densities, respectively.
%
% Examples
% --------
% results = fpk_omega_probability_modal_map([],5e-3);
% results = fpk_omega_probability_modal_map(0.135:5e-4:0.165,5e-3);
%
% IMPORTANT
% ---------
% QUICK mode uses the same weak reset regularization as the original code.
% It is suitable for qualitative identification of peak number and peak
% locations, but not for final relative modal weights or occupation
% probabilities. For publication, verify representative points close to
% every detected transition by strict FPK and/or original-system MC.

clc;
close all;
%save_dir = '不同Omega的概率密度图'; 
%% 1. INPUTS
if nargin < 1 || isempty(Omega_scan)
    % 61 points. Reduce the step to 2.5e-4 after the preliminary scan if a
    % smoother transition boundary is needed.
    Omega_scan = 0.135:5.0e-4:0.165;
end
if nargin < 2 || isempty(S_eta)
    S_eta = 5.0e-3;
end

validateattributes(Omega_scan,{'numeric'}, ...
    {'real','finite','vector','positive'},mfilename,'Omega_scan',1);
validateattributes(S_eta,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'S_eta',2);
Omega_scan = Omega_scan(:).';

%% 2. TLS SLOW-FLOW PARAMETERS
par.k1 = 0.0313;
par.k2 = -0.0609;
par.k3 = 0.0381;
par.b0 = 0.5000;
par.b2 = 0.5417;
par.xi1 = 0.005;
par.xi2 = 0.005;
par.epsv = 0.1;
par.Zhat = 0.001;
par.ghat = 0.001;

par.wn = sqrt(par.k1);
par.rhohat = par.xi1 + par.xi2*par.b0;
par.rhos = par.rhohat/(par.epsv^2);
par.Zs = par.Zhat/(par.epsv^3);
par.gs = par.ghat/(par.epsv^3);
par.kappa = ...
    (3*par.k3-10*par.k2^2/(3*par.wn^2))/(2*par.wn);

%% 3. NUMERICAL AND PEAK-DETECTION SETTINGS
% true: qualitative fast map; false: strict stationary solve retained from
% the original program. Strict mode is much slower over a dense Omega scan.
cfg.quick_debug = true;

if cfg.quick_debug
    % A practical grid for a one-parameter scan. After locating the
    % transition intervals, rerun those intervals with a finer grid.
    cfg.Na = 160;
    cfg.Nphi = 720;
    cfg.a_min = 0.02;
    cfg.a_max = 4.20;
else
    cfg.Na = 200;
    cfg.Nphi = 1600;
    cfg.a_min = 0.02;
    cfg.a_max = 4.50;
end

cfg.mass_defect_tol = 1e-10;
cfg.negative_mass_tol = 1e-9;
cfg.residual_tol = 1e-8;
cfg.dual_solution_L1_tol = 2e-5;
cfg.max_abs_Peclet = 1200;
cfg.max_inner_boundary_mass = 5e-5;
cfg.max_outer_boundary_mass = 5e-5;

% Weak-reset regularized solve.
cfg.quick_reset_rate = 1e-3;
cfg.quick_max_reset_rate = 1e-1;
cfg.quick_max_attempts = 3;
cfg.quick_negative_mass_tol = 1e-6;
cfg.seed_sigma_a = 0.16;
cfg.seed_sigma_phi = 0.22;

cfg.gmres_restart = 40;
cfg.gmres_tol = 1e-8;
cfg.gmres_maxit = 12;
cfg.ilu_droptol = 1e-2;
cfg.direct_fallback_maxN = 150000;

% Peak-identification criteria. These thresholds suppress grid-scale
% numerical ripples without merging the well-separated low/high modes.
cfg.peak_smooth_span = 5;
cfg.peak_min_height_ratio = 0.05;
cfg.peak_min_prominence_ratio = 0.02;
cfg.peak_min_separation_a = 0.25;
cfg.max_number_modes = 2;

% Plot settings retained from the original Fig. 8 code.
cfg.font_name = 'Times New Roman';
cfg.font_size = 24;
cfg.line_width = 3;
cfg.axis_line_width = 1.5;

cfg.save_mat = 'Fig8_Omega_probability_modal_map.mat';
cfg.save_png = 'Fig8_Omega_probability_modal_map.png';
cfg.save_pdf = 'Fig8_Omega_probability_modal_map.pdf';

%% 4. GRID
Na = cfg.Na;
Nphi = cfg.Nphi;
da = (cfg.a_max-cfg.a_min)/Na;
dphi = 2*pi/Nphi;
a = cfg.a_min + ((1:Na)-0.5)*da;
phi = ((1:Nphi)-0.5)*dphi;
[A,PHI] = ndgrid(a,phi);

nOmega = numel(Omega_scan);
pa_all = NaN(Na,nOmega);
mode_count_raw = NaN(1,nOmega);
most_probable_a = NaN(1,nOmega);
peak_amplitudes = cell(1,nOmega);
solve_valid = false(1,nOmega);

fprintf('\n====================================================\n');
fprintf('OMEGA-INDUCED STATIONARY PROBABILITY-MODAL MAP\n');
fprintf('====================================================\n');
fprintf('mode      = %s\n',ternary(cfg.quick_debug, ...
    'QUICK REGULARIZED','STRICT DUAL'));
fprintf('Omega     = [%.6f, %.6f], n = %d\n', ...
    Omega_scan(1),Omega_scan(end),nOmega);
fprintf('S_eta     = %.8e\n',S_eta);
fprintf('grid      = Na x Nphi = %d x %d\n',Na,Nphi);
fprintf('a range   = [%.4f, %.4f]\n',cfg.a_min,cfg.a_max);

%% 5. SOLVE THE STATIONARY FPK EQUATION FOR EACH OMEGA
Qslow = par.epsv^2*par.gs^2*S_eta;
V1 = Qslow/(2*par.wn^2);
cell_area = da*dphi;
N = Na*Nphi;
mass_weight = cell_area*ones(N,1);

for k = 1:nOmega
    Omega = Omega_scan(k);
    sigma = (Omega-par.wn)/(par.epsv^2);

    fprintf('[%3d/%3d] Omega = %.8f ... ',k,nOmega,Omega);

    W1 = ...
        -par.rhos.*A ...
        -(par.xi2*par.b2/4).*A.^3 ...
        +Qslow./(4*A*par.wn^2) ...
        +par.Zs/(2*par.wn).*sin(PHI);

    W2 = ...
        sigma ...
        -par.kappa.*A.^2/4 ...
        +par.Zs./(2*par.wn*A).*cos(PHI);

    V2 = Qslow./(2*par.wn^2*A.^2);

    try
        [LFV,fv_info] = build_FV_matrix_sg( ...
            Na,Nphi,W1,W2,V1,V2,da,dphi);

        mass_defect = norm(mass_weight.'*LFV,inf);
        if mass_defect > cfg.mass_defect_tol
            warning('Mass defect %.3e at Omega=%.8f.',mass_defect,Omega);
            fprintf('FAILED (mass conservation)\n');
            continue;
        end

        if fv_info.max_abs_Peclet > cfg.max_abs_Peclet
            warning('Large max |Pe|=%.2f at Omega=%.8f.', ...
                fv_info.max_abs_Peclet,Omega);
        end

        det_info = deterministic_fixed_points(Omega,par);

        if cfg.quick_debug
            q = build_reference_density(A,PHI,det_info,cell_area,cfg);
            [p,solve_info] = solve_stationary_regularized( ...
                LFV,q,cell_area,cfg);
        else
            [p,solve_info] = solve_stationary_dual( ...
                LFV,Na,Nphi,cell_area,cfg);
        end

        if ~solve_info.valid
            fprintf('FAILED (stationary solve)\n');
            continue;
        end

        P = reshape(p,[Na,Nphi]);
        pa = sum(P,2)*dphi;
        pa = pa/(sum(pa)*da);

        inner_mass = sum(pa(1:min(2,Na)))*da;
        outer_mass = sum(pa(max(1,Na-1):Na))*da;
        if inner_mass > cfg.max_inner_boundary_mass || ...
                outer_mass > cfg.max_outer_boundary_mass
            warning(['Amplitude-domain tail mass is non-negligible at ', ...
                'Omega=%.8f: inner=%.3e, outer=%.3e.'], ...
                Omega,inner_mass,outer_mass);
        end

        [peak_idx,peak_a] = detect_probability_modes(pa,a,da,cfg);
        [~,idx_global] = max(pa);

        pa_all(:,k) = pa;
        mode_count_raw(k) = numel(peak_idx);
        most_probable_a(k) = a(idx_global);
        peak_amplitudes{k} = peak_a(:).';
        solve_valid(k) = true;

        fprintf('modes=%d, a_peak=',numel(peak_idx));
        fprintf(' %.4f',peak_a);
        fprintf(', max|Pe|=%.1f\n',fv_info.max_abs_Peclet);

    catch ME
        warning('Omega=%.8f failed: %s',Omega,ME.message);
        fprintf('FAILED (exception)\n');
    end
end

%% 6. REMOVE ISOLATED ONE-GRID CLASSIFICATION FLIPS
mode_count = clean_isolated_mode_flips(mode_count_raw,solve_valid);
transition_Omega = locate_modal_transitions(Omega_scan,mode_count,solve_valid);

% Reconstruct peak arrays after the mode-count cleaning. If a point was
% reclassified as unimodal, retain only the global mode. If it was
% reclassified as bimodal, retain the two stored peaks when available.
for k = 1:nOmega
    if ~solve_valid(k)
        continue;
    end
    if mode_count(k)==1
        peak_amplitudes{k} = most_probable_a(k);
    elseif mode_count(k)==2 && numel(peak_amplitudes{k})>2
        peak_amplitudes{k} = peak_amplitudes{k}(1:2);
    end
end

%% 7. COLUMN-NORMALIZED PROBABILITY MAP
% Column normalization emphasizes modal topology. It must not be used to
% compare absolute density magnitudes at different Omega values.
pa_normalized = pa_all;
for k = 1:nOmega
    if solve_valid(k)
        vmax = max(pa_all(:,k));
        if isfinite(vmax) && vmax>0
            pa_normalized(:,k) = pa_all(:,k)/vmax;
        end
    end
end

%% 8. FIGURE: OMEGA-INDUCED PROBABILITY-MODAL STRUCTURE
% fig = figure('Units','centimeters','Position',[2 4 36 15]);
% tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

% ----------------------------------------------------------------------
% (a) Stationary normalized amplitude-density map
% ----------------------------------------------------------------------
%% =====================================================================
% Figure 1: Omega-induced stationary probability-modal density map
% =====================================================================

fig1 = figure( ...
    'Units','centimeters', ...
    'Position',[5 5 20 15]);

% 必须先建立坐标轴
ax1 = axes('Parent',fig1);

% pa_normalized的尺寸应为 Na × nOmega
imagesc(ax1,Omega_scan,a,pa_normalized);

% imagesc默认反转纵坐标，这里恢复正常方向
set(ax1,'YDir','normal');

hold(ax1,'on');

% 单独指定当前坐标轴的colormap
colormap(ax1,parula(256));

% 若MATLAB版本较旧，使用caxis代替clim
caxis(ax1,[0,1]);

% Overlay every detected probability-mode location
for k = 1:nOmega

    if ~solve_valid(k) || isempty(peak_amplitudes{k})
        continue;
    end

    current_peaks = peak_amplitudes{k};

    plot(ax1, ...
        Omega_scan(k)*ones(size(current_peaks)), ...
        current_peaks, ...
        'k.', ...
        'MarkerSize',6, ...
        'HandleVisibility','off');
end

% P-bifurcation boundaries: mode number changes 1 <-> 2
for j = 1:numel(transition_Omega)

    xline(ax1,transition_Omega(j), ...
        'k--', ...
        'LineWidth',2.0, ...
        'HandleVisibility','off');
end

xlabel(ax1,'$\Omega$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.font_size, ...
    'FontName',cfg.font_name);

ylabel(ax1,'$a$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.font_size, ...
    'FontName',cfg.font_name);

% Colorbar
cb1 = colorbar(ax1);

ylabel(cb1,'$p(a;\Omega)/\max_a p(a;\Omega)$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.font_size, ...
    'FontName',cfg.font_name);

set(cb1, ...
    'TickLabelInterpreter','latex', ...
    'FontName',cfg.font_name, ...
    'FontSize',cfg.font_size, ...
    'LineWidth',cfg.axis_line_width);

set(ax1, ...
    'FontName',cfg.font_name, ...
    'FontSize',cfg.font_size, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',cfg.axis_line_width, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');

xlim(ax1,[Omega_scan(1),Omega_scan(end)]);
ylim(ax1,[cfg.a_min,cfg.a_max]);

grid(ax1,'on');
box(ax1,'on');
save_all_formats(fig1, save_dir, '6_Omega_pa_pphi_map');
% 可选：子图编号
% text(ax1,0.5,-0.20,'(a)', ...
%     'Units','normalized', ...
%     'HorizontalAlignment','center', ...
%     'Interpreter','latex', ...
%     'FontName',cfg.font_name, ...
%     'FontSize',cfg.font_size);


%% =====================================================================
% Figure 2: Detected probability-modal amplitudes
% =====================================================================

fig2 = figure( ...
    'Units','centimeters', ...
    'Position',[5 5 20 15]);

% 必须建立第二个坐标轴
ax2 = axes('Parent',fig2);

hold(ax2,'on');

% 用于图例的虚拟对象
h_uni = plot(ax2,NaN,NaN, ...
    'o', ...
    'Color',[0 0 1], ...
    'MarkerFaceColor',[0 0 1], ...
    'MarkerSize',7, ...
    'LineWidth',1.2);

h_bi = plot(ax2,NaN,NaN, ...
    'o', ...
    'Color',[1 0 0], ...
    'MarkerFaceColor',[1 0 0], ...
    'MarkerSize',7, ...
    'LineWidth',1.2);

% Plot detected modal amplitudes
for k = 1:nOmega

    if ~solve_valid(k) || isempty(peak_amplitudes{k})
        continue;
    end

    current_peaks = peak_amplitudes{k};

    if mode_count(k) == 1
        marker_color = [0 0 1];
    elseif mode_count(k) >= 2
        marker_color = [1 0 0];
    else
        continue;
    end

    plot(ax2, ...
        Omega_scan(k)*ones(size(current_peaks)), ...
        current_peaks, ...
        'o', ...
        'Color',marker_color, ...
        'MarkerFaceColor',marker_color, ...
        'MarkerSize',6, ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');
end

% P-bifurcation boundaries
for j = 1:numel(transition_Omega)

    xline(ax2,transition_Omega(j), ...
        'k--', ...
        'LineWidth',2.0, ...
        'HandleVisibility','off');
end

xlabel(ax2,'$\Omega$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.font_size, ...
    'FontName',cfg.font_name);

ylabel(ax2,'$a_{\mathrm{peak}}$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.font_size, ...
    'FontName',cfg.font_name);

legend(ax2,[h_uni,h_bi], ...
    {'Unimodal','Bimodal'}, ...
    'Location','best', ...
    'FontSize',cfg.font_size, ...
    'FontName',cfg.font_name, ...
    'Interpreter','latex', ...
    'Box','off');

set(ax2, ...
    'FontName',cfg.font_name, ...
    'FontSize',cfg.font_size, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',cfg.axis_line_width, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');

xlim(ax2,[Omega_scan(1),Omega_scan(end)]);
ylim(ax2,[cfg.a_min,cfg.a_max]);

grid(ax2,'on');
box(ax2,'on');
save_all_formats(fig2, save_dir, '6_Omega_mode');
% 可选：子图编号
% text(ax2,0.5,-0.20,'(b)', ...
%     'Units','normalized', ...
%     'HorizontalAlignment','center', ...
%     'Interpreter','latex', ...
%     'FontName',cfg.font_name, ...
%     'FontSize',cfg.font_size);

% %% 9. SAVE RESULTS
% results = struct();
% results.parameters = par;
% results.settings = cfg;
% results.Omega = Omega_scan(:);
% results.S_eta = S_eta;
% results.a = a(:);
% results.pa = pa_all;
% results.pa_normalized = pa_normalized;
% results.solve_valid = solve_valid(:);
% results.mode_count_raw = mode_count_raw(:);
% results.mode_count = mode_count(:);
% results.most_probable_a = most_probable_a(:);
% results.peak_amplitudes = peak_amplitudes;
% results.transition_Omega = transition_Omega(:);
% 
% save(cfg.save_mat,'results','-v7.3');
% try
%     exportgraphics(fig,cfg.save_png,'Resolution',600);
%     exportgraphics(fig,cfg.save_pdf,'ContentType','vector');
% catch
%     print(fig,cfg.save_png,'-dpng','-r600');
%     print(fig,cfg.save_pdf,'-dpdf','-painters');
% end
% 
% fprintf('\nDetected P-bifurcation boundaries at Omega =');
% if isempty(transition_Omega)
%     fprintf(' none');
% else
%     fprintf(' %.8f',transition_Omega);
% end
% fprintf('\nSaved: %s\nSaved: %s\nSaved: %s\n', ...
%     cfg.save_mat,cfg.save_png,cfg.save_pdf);

end

%% =====================================================================
% Significant probability-mode detection
% =====================================================================
function [peak_idx,peak_a] = detect_probability_modes(pa,a,da,cfg)

y = pa(:);
y(~isfinite(y)) = 0;
y(y<0) = 0;

span = max(1,round(cfg.peak_smooth_span));
if mod(span,2)==0
    span = span+1;
end
if span>1
    kernel = ones(span,1)/span;
    ys = conv(y,kernel,'same');
else
    ys = y;
end

if max(ys)<=0
    peak_idx = [];
    peak_a = [];
    return;
end

candidate = find(ys(2:end-1)>ys(1:end-2) & ...
                 ys(2:end-1)>=ys(3:end))+1;

% Always retain the global maximum as a candidate.
[~,idx_global] = max(ys);
candidate = unique([candidate(:);idx_global]);

height_threshold = cfg.peak_min_height_ratio*max(ys);
half_window = max(2,round(cfg.peak_min_separation_a/da));
keep = false(size(candidate));
prominence = zeros(size(candidate));

for r = 1:numel(candidate)
    idx = candidate(r);
    iL = max(1,idx-half_window);
    iR = min(numel(ys),idx+half_window);
    left_min = min(ys(iL:idx));
    right_min = min(ys(idx:iR));
    prominence(r) = ys(idx)-max(left_min,right_min);
    keep(r) = ys(idx)>=height_threshold && ...
        prominence(r)>=cfg.peak_min_prominence_ratio*max(ys);
end

candidate = candidate(keep);
if isempty(candidate)
    candidate = idx_global;
end

% Enforce a minimum amplitude separation by accepting higher peaks first.
[~,order] = sort(ys(candidate),'descend');
selected = [];
min_sep_idx = max(1,round(cfg.peak_min_separation_a/da));
for r = order(:).'
    idx = candidate(r);
    if isempty(selected) || all(abs(idx-selected)>=min_sep_idx)
        selected(end+1) = idx; %#ok<AGROW>
    end
    if numel(selected)>=cfg.max_number_modes
        break;
    end
end

peak_idx = sort(selected(:));
peak_a = a(peak_idx);
end

%% =====================================================================
% Remove isolated one-grid flips in mode classification
% =====================================================================
function mode_count = clean_isolated_mode_flips(mode_count_raw,valid)
mode_count = mode_count_raw;
for pass = 1:2
    old = mode_count;
    for k = 2:numel(mode_count)-1
        if valid(k-1) && valid(k) && valid(k+1) && ...
                isfinite(old(k-1)) && isfinite(old(k)) && ...
                isfinite(old(k+1)) && ...
                old(k-1)==old(k+1) && old(k)~=old(k-1)
            mode_count(k) = old(k-1);
        end
    end
end
end

%% =====================================================================
% Locate 1 <-> 2 modal-transition boundaries
% =====================================================================
function transition_Omega = locate_modal_transitions(Omega,mode_count,valid)
transition_Omega = [];
for k = 1:numel(Omega)-1
    if valid(k) && valid(k+1) && ...
            isfinite(mode_count(k)) && isfinite(mode_count(k+1)) && ...
            mode_count(k)~=mode_count(k+1)
        transition_Omega(end+1) = 0.5*(Omega(k)+Omega(k+1)); %#ok<AGROW>
    end
end
end

%% =====================================================================
% Conservative Scharfetter-Gummel finite-volume matrix
% =====================================================================
function [LFV,info] = build_FV_matrix_sg( ...
    Na,Nphi,W1,W2,V1,V2,da,dphi)

N = Na*Nphi;
number_a_faces = (Na-1)*Nphi;
number_phi_faces = Na*Nphi;
max_entries = 4*(number_a_faces+number_phi_faces);

I = zeros(max_entries,1);
J = zeros(max_entries,1);
S = zeros(max_entries,1);
counter = 0;
index = @(i,j) i+(j-1)*Na;

Wface_a = 0.5*(W1(1:Na-1,:)+W1(2:Na,:));
D_a = 0.5*V1;
Pe_a = Wface_a*da/D_a;
coeffL_a = (D_a/da).*bernoulli_sg(-Pe_a);
coeffR_a = -(D_a/da).*bernoulli_sg(Pe_a);

for j = 1:Nphi
    for i = 1:Na-1
        cellL = index(i,j);
        cellR = index(i+1,j);
        coeffL = coeffL_a(i,j);
        coeffR = coeffR_a(i,j);

        counter=counter+1; I(counter)=cellL; J(counter)=cellL;
        S(counter)=-coeffL/da;
        counter=counter+1; I(counter)=cellL; J(counter)=cellR;
        S(counter)=-coeffR/da;
        counter=counter+1; I(counter)=cellR; J(counter)=cellL;
        S(counter)= coeffL/da;
        counter=counter+1; I(counter)=cellR; J(counter)=cellR;
        S(counter)= coeffR/da;
    end
end

jp_all = [2:Nphi,1];
Wface_phi = 0.5*(W2+W2(:,jp_all));
D_phi = 0.25*(V2+V2(:,jp_all));
Pe_phi = Wface_phi*dphi./D_phi;
coeffL_phi = (D_phi/dphi).*bernoulli_sg(-Pe_phi);
coeffR_phi = -(D_phi/dphi).*bernoulli_sg(Pe_phi);

for j = 1:Nphi
    jp = j+1;
    if jp > Nphi
        jp = 1;
    end

    for i = 1:Na
        cellL = index(i,j);
        cellR = index(i,jp);
        coeffL = coeffL_phi(i,j);
        coeffR = coeffR_phi(i,j);

        counter=counter+1; I(counter)=cellL; J(counter)=cellL;
        S(counter)=-coeffL/dphi;
        counter=counter+1; I(counter)=cellL; J(counter)=cellR;
        S(counter)=-coeffR/dphi;
        counter=counter+1; I(counter)=cellR; J(counter)=cellL;
        S(counter)= coeffL/dphi;
        counter=counter+1; I(counter)=cellR; J(counter)=cellR;
        S(counter)= coeffR/dphi;
    end
end

LFV = sparse(I(1:counter),J(1:counter),S(1:counter),N,N);

info = struct();
info.max_abs_Peclet = max([abs(Pe_a(:));abs(Pe_phi(:))]);
info.max_abs_Peclet_a = max(abs(Pe_a(:)));
info.max_abs_Peclet_phi = max(abs(Pe_phi(:)));
end

function B = bernoulli_sg(z)
% Stable B(z)=z/(exp(z)-1).
B = zeros(size(z));
small = abs(z)<1e-6;
large_pos = z>50;
large_neg = z<-50;
mid = ~(small | large_pos | large_neg);

zs = z(small);
B(small) = 1-zs/2+zs.^2/12-zs.^4/720;
B(large_pos) = z(large_pos).*exp(-z(large_pos));
B(large_neg) = -z(large_neg);
B(mid) = z(mid)./expm1(z(mid));
end

%% =====================================================================
% Quick unique stationary solve with weak reset regularization
% =====================================================================
function [p,info] = solve_stationary_regularized( ...
    LFV,q,cell_area,cfg)

N = size(LFV,1);
q = q(:);
q(q<0) = 0;
q = q/(sum(q)*cell_area);

info = initialize_solve_info();
info.solver_mode = 'quick-reset-GMRES';

alpha = cfg.quick_reset_rate;
last_message = '';

for attempt = 1:cfg.quick_max_attempts
    Areg = alpha*speye(N)-LFV;
    rhs = alpha*q;
    x0 = q;

    [p_try,flag,relres,iter,method_name,last_message] = ...
        solve_shifted_linear_system(Areg,rhs,x0,cfg);

    if all(isfinite(p_try))
        negative_mass = max(-sum(min(p_try,0))*cell_area,0);
        if negative_mass <= cfg.quick_negative_mass_tol
            p_try(p_try<0) = 0;
            mass = sum(p_try)*cell_area;
            if isfinite(mass) && mass>0
                p_try = p_try/mass;

                reg_residual = norm(Areg*p_try-rhs,inf);
                fpk_residual = norm(LFV*p_try,inf);

                % Accept either a formally converged iterative solve or a
                % very small explicit residual after normalization.
                accepted = (flag==0) || ...
                    (relres<=10*cfg.gmres_tol) || ...
                    (reg_residual<=1e-8*max(1,norm(rhs,inf)));

                if accepted
                    p = p_try;
                    info.valid = true;
                    info.negative_mass = negative_mass;
                    info.residual_Linf = fpk_residual;
                    info.regularized_residual_Linf = reg_residual;
                    info.linear_flag = flag;
                    info.linear_relres = relres;
                    info.linear_iterations = iter_to_scalar(iter,cfg.gmres_restart);
                    info.reset_rate = alpha;
                    info.solver_mode = method_name;
                    return;
                end
            end
        end
    end

    alpha = min(10*alpha,cfg.quick_max_reset_rate);
end

p = NaN(N,1);
warning('Regularized solve failed. Last linear-solver message: %s', ...
    last_message);
end

function [x,flag,relres,iter,method_name,message] = ...
    solve_shifted_linear_system(A,b,x0,cfg)

N = size(A,1);
x = NaN(N,1);
flag = Inf;
relres = Inf;
iter = [0,0];
method_name = 'none';
message = '';

% First choice: ILUTP-preconditioned GMRES.
try
    setup.type = 'ilutp';
    setup.droptol = cfg.ilu_droptol;
    setup.udiag = 1;
    [L,U] = ilu(A,setup);

    [x,flag,relres,iter] = gmres(A,b,cfg.gmres_restart, ...
        cfg.gmres_tol,cfg.gmres_maxit,L,U,x0);
    method_name = 'quick-reset-GMRES-ILUTP';

    if all(isfinite(x)) && ...
            (flag==0 || relres<=10*cfg.gmres_tol)
        return;
    end
catch ME
    message = ME.message;
end

% Second choice: Jacobi-preconditioned GMRES, cheaper and robust to ILU
% breakdowns caused by extremely anisotropic transition rates.
try
    d = diag(A);
    tiny = abs(d)<eps(max(1,max(abs(d))));
    d(tiny) = 1;
    M = spdiags(d,0,N,N);

    [x,flag,relres,iter] = gmres(A,b,cfg.gmres_restart, ...
        cfg.gmres_tol,2*cfg.gmres_maxit,M,[],x0);
    method_name = 'quick-reset-GMRES-Jacobi';

    if all(isfinite(x)) && ...
            (flag==0 || relres<=10*cfg.gmres_tol)
        return;
    end
catch ME
    message = ME.message;
end

% Final fallback.  The reset shift makes A nonsingular, unlike the original
% row-replaced null-space solve.  This branch is used only if iterative
% methods fail and the matrix is not too large.
if N <= cfg.direct_fallback_maxN
    try
        x = A\b;
        flag = 0;
        relres = norm(A*x-b)/max(norm(b),realmin);
        iter = [0,0];
        method_name = 'quick-reset-direct';
        return;
    catch ME
        message = ME.message;
    end
end
end

%% =====================================================================
% Reference density for quick mode
% =====================================================================
function q = build_reference_density(A,PHI,det_info,cell_area,cfg)

stable_idx = find(det_info.stable);

if isempty(stable_idx)
    % Broad positive reference if no stable deterministic root is detected.
    q = ones(size(A));
else
    q = zeros(size(A));
    for r = stable_idx(:).'
        ar = det_info.a_roots(r);
        phir = det_info.phi_roots(r);
        dph = atan2(sin(PHI-phir),cos(PHI-phir));

        qr = exp(-0.5*((A-ar)/cfg.seed_sigma_a).^2 ...
            -0.5*(dph/cfg.seed_sigma_phi).^2);
        q = q + qr;
    end
end

q = max(q,realmin);
q = q/(sum(q(:))*cell_area);
end

%% =====================================================================
% Strict dual-normalization solver retained for optional final calculations
% =====================================================================
function [p,info] = solve_stationary_dual( ...
    LFV,Na,Nphi,cell_area,cfg)

N = size(LFV,1);
diagL = abs(full(diag(LFV)));
[~,row1] = max(diagL);

i2 = max(1,min(Na,round(0.55*Na)));
j2 = max(1,min(Nphi,round(0.37*Nphi)));
row2 = i2+(j2-1)*Na;
if row2==row1
    row2 = max(1,min(N,row1+round(N/3)));
end

[p1,info1] = solve_one_normalization_row( ...
    LFV,row1,cell_area,cfg);
[p2,info2] = solve_one_normalization_row( ...
    LFV,row2,cell_area,cfg);

info = initialize_solve_info();
info.solver_mode = 'strict-dual-direct';
info.row1 = row1;
info.row2 = row2;

if ~(info1.valid && info2.valid)
    p = NaN(N,1);
    return;
end

dualL1 = sum(abs(p1-p2))*cell_area;
p = 0.5*(p1+p2);
p = max(p,0);
p = p/(sum(p)*cell_area);

negative_mass = max(info1.negative_mass,info2.negative_mass);
residual_Linf = max(abs(LFV*p));

info.valid = ...
    dualL1<=cfg.dual_solution_L1_tol && ...
    negative_mass<=cfg.negative_mass_tol && ...
    residual_Linf<=cfg.residual_tol && ...
    all(isfinite(p));
info.negative_mass = negative_mass;
info.residual_Linf = residual_Linf;
info.dual_solution_L1 = dualL1;
info.linear_flag = 0;
info.linear_relres = residual_Linf;
end

function [p,info] = solve_one_normalization_row( ...
    LFV,replace_row,cell_area,cfg)

N = size(LFV,1);
Aeq = LFV;
rhs = zeros(N,1);
Aeq(replace_row,:) = sparse(1,1:N,ones(1,N),1,N);
rhs(replace_row) = 1/cell_area;

try
    p = Aeq\rhs;
catch
    p = NaN(N,1);
end

info = initialize_solve_info();

if any(~isfinite(p))
    return;
end

negative_mass = max(-sum(min(p,0))*cell_area,0);
if negative_mass > cfg.negative_mass_tol
    return;
end
p(p<0) = 0;

mass = sum(p)*cell_area;
if ~isfinite(mass) || mass<=0
    return;
end
p = p/mass;

residual_Linf = max(abs(LFV*p));
info.valid = ...
    negative_mass<=cfg.negative_mass_tol && ...
    residual_Linf<=cfg.residual_tol;
info.negative_mass = negative_mass;
info.residual_Linf = residual_Linf;
end

%% =====================================================================
% Deterministic slow-flow fixed points
% =====================================================================
function info = deterministic_fixed_points(Omega,par)

sigma = (Omega-par.wn)/(par.epsv^2);
c = par.xi2*par.b2/4;
d = par.kappa/4;
F = par.Zs/(2*par.wn);

coef = [ ...
    c^2+d^2, ...
    2*(par.rhos*c-sigma*d), ...
    par.rhos^2+sigma^2, ...
    -F^2];

yroot = roots(coef);
mask = abs(imag(yroot))<1e-8 & real(yroot)>0;
a_roots = sort(sqrt(real(yroot(mask))));

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

    eigJ = eig([J11 J12;J21 J22]);
    max_real_eig(j) = max(real(eigJ));
    stable(j) = max_real_eig(j)<0;
    phi_roots(j) = phi_j;
end

info = struct();
info.a_roots = a_roots;
info.phi_roots = phi_roots;
info.stable = stable;
info.max_real_eig = max_real_eig;
end

%% =====================================================================
% Output helpers
% =====================================================================
function write_marginal_tables(a,phi,pa,pphi,Omega_list,pa_file,pphi_file)

var_names = cell(1,numel(Omega_list)+1);
var_names{1} = 'a';
for k = 1:numel(Omega_list)
    var_names{k+1} = matlab.lang.makeValidName( ...
        sprintf('p_a_Omega_%.6f',Omega_list(k)));
end
Tpa = array2table([a(:),pa],'VariableNames',var_names);
writetable(Tpa,pa_file);

var_names_phi = cell(1,numel(Omega_list)+1);
var_names_phi{1} = 'phi';
for k = 1:numel(Omega_list)
    var_names_phi{k+1} = matlab.lang.makeValidName( ...
        sprintf('p_phi_Omega_%.6f',Omega_list(k)));
end
Tphi = array2table([phi(:),pphi],'VariableNames',var_names_phi);
writetable(Tphi,pphi_file);
end

function s = latex_scientific(x)
if x==0
    s = '0';
    return;
end
exponent = floor(log10(abs(x)));
mantissa = x/10^exponent;
if abs(mantissa-round(mantissa))<1e-12
    mantissa_text = sprintf('%d',round(mantissa));
else
    mantissa_text = sprintf('%.3g',mantissa);
end
s = sprintf('%s\\times10^{%d}',mantissa_text,exponent);
end

function d = empty_diagnostics()
d = struct( ...
    'Omega',NaN, ...
    'S_eta',NaN, ...
    'Ea',NaN, ...
    'Ea2',NaN, ...
    'Var_a',NaN, ...
    'mode_a',NaN, ...
    'mode_phi',NaN, ...
    'joint_mass',NaN, ...
    'pa_mass',NaN, ...
    'pphi_mass',NaN, ...
    'mass_defect',NaN, ...
    'negative_mass',NaN, ...
    'residual_Linf',NaN, ...
    'dual_solution_L1',NaN, ...
    'max_abs_Peclet',NaN, ...
    'max_abs_Peclet_a',NaN, ...
    'max_abs_Peclet_phi',NaN, ...
    'inner_boundary_mass',NaN, ...
    'outer_boundary_mass',NaN, ...
    'deterministic_amplitudes',[], ...
    'deterministic_phases',[], ...
    'deterministic_stability',[], ...
    'solver_mode','', ...
    'reset_rate',NaN, ...
    'linear_flag',NaN, ...
    'linear_relres',NaN, ...
    'linear_iterations',NaN);
end

function info = initialize_solve_info()
info = struct();
info.valid = false;
info.negative_mass = Inf;
info.residual_Linf = Inf;
info.regularized_residual_Linf = Inf;
info.dual_solution_L1 = NaN;
info.row1 = NaN;
info.row2 = NaN;
info.solver_mode = '';
info.reset_rate = NaN;
info.linear_flag = Inf;
info.linear_relres = Inf;
info.linear_iterations = NaN;
end

function n = iter_to_scalar(iter,restart)
if isempty(iter)
    n = NaN;
elseif isscalar(iter)
    n = iter;
else
    n = iter(1)*restart + iter(2);
end
end

function out = ternary(condition,a,b)
if condition
    out = a;
else
    out = b;
end
end
