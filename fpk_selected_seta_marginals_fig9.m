function results = fpk_selected_seta_marginals_fig9(S_eta_list)
%FPK_SELECTED_SETA_MARGINALS_FIG9
% Fast/robust stationary FPK solver for qualitative inspection of the
% single-peak/double-peak evolution of p(a) and p(phi).
%
% QUICK MODE (default)
% --------------------
% The low-noise bistable FPK matrix can be numerically almost reducible.
% Replacing one matrix row by the normalization condition may then return a
% singular/ill-conditioned solve or two mutually disjoint null vectors.
% Quick mode solves the weakly regularized stationary problem
%
%     0 = L*p + alpha*(q-p),
%
% where q is an equal-weight reference density centered at the stable
% deterministic slow-flow fixed points.  The small reset rate alpha makes
% the linear problem unique and is intended only for rapid qualitative
% identification of peak number and peak locations.  Relative peak weights,
% means and variances from quick mode must not be used as final publication
% values.
%
% Examples
% --------
% S_eta_list = [0.0001,0.0005,0.001,0.005,0.01];
% results = fpk_selected_seta_marginals_fig9(S_eta_list);
%
% Omega=0.1530 and xi1=xi2=0.005 are fixed.

clc;
close all;

save_dir = '不同S_eta的概率密度图';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% 1. INPUTS
if nargin < 1 || isempty(S_eta_list)
    S_eta_list = [0.001,0.003,0.005,0.007,0.009];
end

validateattributes(S_eta_list,{'numeric'}, ...
    {'real','finite','vector','positive'},mfilename,'S_eta_list',1);

S_eta_list = S_eta_list(:).';
Omega = 0.1530;

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

%% 3. SETTINGS
% Keep true while inspecting the noise-intensity-dependent peak topology.
cfg.quick_debug = true;

if cfg.quick_debug
    % Grid used for the five selected noise intensities.
    cfg.Na = 200;
    cfg.Nphi = 2400;
    cfg.a_min = 0.02;
    cfg.a_max = 4.70;
else
    % Strict mode retains the original dual-normalization solve.  It may be
    % ill-conditioned in the very-low-noise bistable regime.
    cfg.Na = 200;
    cfg.Nphi = 1600;
    cfg.a_min = 0.02;
    cfg.a_max = 4.70;
end

cfg.mass_defect_tol = 1e-10;
cfg.negative_mass_tol = 1e-9;
cfg.residual_tol = 1e-8;
cfg.dual_solution_L1_tol = 2e-5;
cfg.max_abs_Peclet = 600;

cfg.max_inner_boundary_mass = 5e-5;
cfg.max_outer_boundary_mass = 5e-5;

% Fast regularized solve settings.
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

cfg.line_width = 3;
cfg.save_mat = 'Fig9_selected_Seta_stationary_FPK_fast.mat';
cfg.save_pa_csv = 'Fig9_selected_Seta_pa_fast.csv';
cfg.save_pphi_csv = 'Fig9_selected_Seta_pphi_fast.csv';
cfg.save_moments_csv = 'Fig9_selected_Seta_moments_fast.csv';
cfg.save_figure_pa = 'Fig9_Seta_induced_pa_fast.png';
cfg.save_figure_phi = 'Fig9_Seta_induced_pphi_fast.png';
cfg.save_figure_pa_pdf = 'Fig9_Seta_induced_pa_fast.pdf';
cfg.save_figure_phi_pdf = 'Fig9_Seta_induced_pphi_fast.pdf';

%% 4. GRID
Na = cfg.Na;
Nphi = cfg.Nphi;
da = (cfg.a_max-cfg.a_min)/Na;
dphi = 2*pi/Nphi;
a = cfg.a_min + ((1:Na)-0.5)*da;
phi = ((1:Nphi)-0.5)*dphi;
[A,PHI] = ndgrid(a,phi);

nSeta = numel(S_eta_list);
pa_all = NaN(Na,nSeta);
pphi_all = NaN(Nphi,nSeta);
joint_pdf = cell(1,nSeta);
diagnostics = repmat(empty_diagnostics(),1,nSeta);

fprintf('\n====================================================\n');
fprintf('STATIONARY FPK AT SELECTED NOISE INTENSITIES\n');
fprintf('====================================================\n');
fprintf('mode   = %s\n',ternary(cfg.quick_debug,'QUICK REGULARIZED','STRICT DUAL'));
fprintf('Omega  = %.8f\n',Omega);
fprintf('xi1    = %.8f\n',par.xi1);
fprintf('xi2    = %.8f\n',par.xi2);
fprintf('S_eta values =');
fprintf(' %.4g',S_eta_list);
fprintf('\n');
fprintf('g_hat  = %.8e\n',par.ghat);
fprintf('grid   = Na x Nphi = %d x %d\n',Na,Nphi);
fprintf('a range= [%.4f, %.4f]\n',cfg.a_min,cfg.a_max);

%% 5. SOLVE EACH NOISE INTENSITY
sigma = (Omega-par.wn)/(par.epsv^2);

% Deterministic fixed points do not depend on S_eta and are unchanged
% throughout this sweep.
det_info_fixed = deterministic_fixed_points(Omega,par);

for k = 1:nSeta
    S_eta = S_eta_list(k);

    fprintf('\n[%d/%d] S_eta = %.8e, Omega = %.6f\n', ...
        k,nSeta,S_eta,Omega);

    % The eps^2 factors are retained because S_eta is the same original-time
    % noise-intensity parameter in the original MC system and in the FPK code.
    Qslow = par.epsv^2*par.gs^2*S_eta;

    W1 = ...
        -par.rhos.*A ...
        -(par.xi2*par.b2/4).*A.^3 ...
        +Qslow./(4*A*par.wn^2) ...
        +par.Zs/(2*par.wn).*sin(PHI);

    W2 = ...
        sigma ...
        -par.kappa.*A.^2/4 ...
        +par.Zs./(2*par.wn*A).*cos(PHI);

    V1 = Qslow/(2*par.wn^2);
    V2 = Qslow./(2*par.wn^2*A.^2);

    [LFV,fv_info] = build_FV_matrix_sg( ...
        Na,Nphi,W1,W2,V1,V2,da,dphi);

    N = Na*Nphi;
    cell_area = da*dphi;
    mass_weight = cell_area*ones(N,1);
    mass_defect = norm(mass_weight.'*LFV,inf);

    fprintf('  pre-solve max |Pe_a|   = %.6f\n', ...
        fv_info.max_abs_Peclet_a);
    fprintf('  pre-solve max |Pe_phi| = %.6f\n', ...
        fv_info.max_abs_Peclet_phi);

    if mass_defect > cfg.mass_defect_tol
        error(['Mass conservation failed at Omega=%.10f: ', ...
            'mass defect=%.3e.'],Omega,mass_defect);
    end

    if fv_info.max_abs_Peclet > cfg.max_abs_Peclet
        if cfg.quick_debug
            warning(['Quick grid has max |Pe|=%.2f at S_eta=%.3e. ', ...
                'Continuing for qualitative peak diagnostics; perform ', ...
                'a grid-convergence check before publication.'], ...
                fv_info.max_abs_Peclet,S_eta);
        else
            error(['The grid is numerically unresolved at S_eta=%.3e: ', ...
                'max |Pe|=%.2f exceeds %.2f.'], ...
                S_eta,fv_info.max_abs_Peclet,cfg.max_abs_Peclet);
        end
    end

    det_info = det_info_fixed;

    if cfg.quick_debug
        q = build_reference_density( ...
            A,PHI,det_info,cell_area,cfg);
        [p,solve_info] = solve_stationary_regularized( ...
            LFV,q,cell_area,cfg);
    else
        [p,solve_info] = solve_stationary_dual( ...
            LFV,Na,Nphi,cell_area,cfg);
    end

    if ~solve_info.valid
        error(['Stationary FPK solve failed at Omega=%.10f: ', ...
            'negative mass=%.3e, residual=%.3e, ', ...
            'linear relres=%.3e.'], ...
            Omega,solve_info.negative_mass, ...
            solve_info.residual_Linf,solve_info.linear_relres);
    end

    P = reshape(p,[Na,Nphi]);

    pa = sum(P,2)*dphi;
    pa = pa/(sum(pa)*da);

    pphi = sum(P,1).'*da;
    pphi = pphi/(sum(pphi)*dphi);

    joint_mass = sum(P(:))*cell_area;
    pa_mass = sum(pa)*da;
    pphi_mass = sum(pphi)*dphi;

    inner_boundary_mass = sum(pa(1:min(2,Na)))*da;
    outer_boundary_mass = sum(pa(max(1,Na-1):Na))*da;

    if inner_boundary_mass > cfg.max_inner_boundary_mass || ...
            outer_boundary_mass > cfg.max_outer_boundary_mass
        if cfg.quick_debug
            warning(['Quick-mode boundary warning at Omega=%.10f: ', ...
                'inner=%.3e, outer=%.3e.'], ...
                Omega,inner_boundary_mass,outer_boundary_mass);
        else
            error(['Amplitude domain is too narrow at Omega=%.10f: ', ...
                'inner=%.3e, outer=%.3e.'], ...
                Omega,inner_boundary_mass,outer_boundary_mass);
        end
    end

    Ea = sum(a(:).*pa)*da;
    Ea2 = sum((a(:).^2).*pa)*da;
    Var_a = max(Ea2-Ea^2,0);

    [~,idx_mode_a] = max(pa);
    [~,idx_mode_phi] = max(pphi);
    mode_a = a(idx_mode_a);
    mode_phi = phi(idx_mode_phi);

    pa_all(:,k) = pa;
    pphi_all(:,k) = pphi;
    joint_pdf{k} = P;

    diagnostics(k).Omega = Omega;
    diagnostics(k).xi1 = par.xi1;
    diagnostics(k).xi2 = par.xi2;
    diagnostics(k).S_eta = S_eta;
    diagnostics(k).Ea = Ea;
    diagnostics(k).Ea2 = Ea2;
    diagnostics(k).Var_a = Var_a;
    diagnostics(k).mode_a = mode_a;
    diagnostics(k).mode_phi = mode_phi;
    diagnostics(k).joint_mass = joint_mass;
    diagnostics(k).pa_mass = pa_mass;
    diagnostics(k).pphi_mass = pphi_mass;
    diagnostics(k).mass_defect = mass_defect;
    diagnostics(k).negative_mass = solve_info.negative_mass;
    diagnostics(k).residual_Linf = solve_info.residual_Linf;
    diagnostics(k).dual_solution_L1 = solve_info.dual_solution_L1;
    diagnostics(k).max_abs_Peclet = fv_info.max_abs_Peclet;
    diagnostics(k).max_abs_Peclet_a = fv_info.max_abs_Peclet_a;
    diagnostics(k).max_abs_Peclet_phi = fv_info.max_abs_Peclet_phi;
    diagnostics(k).inner_boundary_mass = inner_boundary_mass;
    diagnostics(k).outer_boundary_mass = outer_boundary_mass;
    diagnostics(k).deterministic_amplitudes = det_info.a_roots;
    diagnostics(k).deterministic_phases = det_info.phi_roots;
    diagnostics(k).deterministic_stability = det_info.stable;
    diagnostics(k).solver_mode = solve_info.solver_mode;
    diagnostics(k).reset_rate = solve_info.reset_rate;
    diagnostics(k).linear_flag = solve_info.linear_flag;
    diagnostics(k).linear_relres = solve_info.linear_relres;
    diagnostics(k).linear_iterations = solve_info.linear_iterations;

    fprintf('  E[a]                    = %.10f\n',Ea);
    fprintf('  Var(a)                  = %.10e\n',Var_a);
    fprintf('  mode of p(a)            = %.10f\n',mode_a);
    fprintf('  mode of p(phi)          = %.10f rad\n',mode_phi);
    fprintf('  integral p(a) da        = %.12f\n',pa_mass);
    fprintf('  integral p(phi) dphi    = %.12f\n',pphi_mass);
    fprintf('  integral p(a,phi)       = %.12f\n',joint_mass);
    fprintf('  max |cell Peclet|       = %.4f\n',fv_info.max_abs_Peclet);
    fprintf('  original FPK residual   = %.3e\n',solve_info.residual_Linf);
    fprintf('  linear solver relres    = %.3e\n',solve_info.linear_relres);
    fprintf('  solver mode             = %s\n',solve_info.solver_mode);
    if cfg.quick_debug
        fprintf('  reset rate alpha        = %.3e\n',solve_info.reset_rate);
    else
        fprintf('  dual-solution L1 diff.  = %.3e\n', ...
            solve_info.dual_solution_L1);
    end
    fprintf('  inner/outer tail mass   = %.3e / %.3e\n', ...
        inner_boundary_mass,outer_boundary_mass);

    if isempty(det_info.a_roots)
        fprintf('  deterministic roots     = none\n');
    else
        fprintf('  deterministic roots a   =');
        fprintf(' %.6f',det_info.a_roots);
        fprintf('\n');
        fprintf('  deterministic stability =');
        fprintf(' %d',det_info.stable);
        fprintf('\n');
    end
end

%% 6. FIG. 9 PLOTS
% Omega=0.1530 and xi1=xi2=0.005 are fixed, while S_eta is varied.

% Deep blue, red, black, cyan, magenta.
color_order = [ ...
       0, 0, 1; ...  % 蓝,
       1, 0, 0; ...  % 红  
       0, 0, 0;...    %黑
       0, 1, 1; ...   %青
       1, 0, 1];     % 紫色 

seta_labels = arrayfun( ...
    @(x) ['$S_{\eta}=',num2str(x,'%.4g'),'$'], ...
    S_eta_list, ...
    'UniformOutput',false);

% ----------------------------------------------------------------------
% Fig. 9(a): p(a)
% ----------------------------------------------------------------------
fig_pa = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 20 15]);

hold on;
h_a = gobjects(nSeta,1);

for k = 1:nSeta
    h_a(k) = plot( ...
        a,pa_all(:,k), ...
        'LineWidth',3, ...
        'Color',color_order(k,:));
end

xlabel('$a$', ...
    'Interpreter','latex', ...
    'FontSize',24, ...
    'FontName','Times New Roman');

ylabel('$p(a)$', ...
    'Interpreter','latex', ...
    'FontSize',24, ...
    'FontName','Times New Roman');

legend( ...
    h_a,seta_labels, ...
    'Location','northwest', ...
    'FontSize',24, ...
    'FontName','Times New Roman', ...
    'Interpreter','latex', ...
    'Box','off');

ax = gca;
set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',24, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.5, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');
xlim([1,4.2]);
%xlim([cfg.a_min,cfg.a_max]);
ylim([0,1.30*max(pa_all(:))]);

grid on;
box on;

save_all_formats(fig_pa,save_dir,'6_guass_Seta_pa');

% ----------------------------------------------------------------------
% Fig. 9(b): p(varphi)
% ----------------------------------------------------------------------
fig_phi = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 20 15]);

hold on;
h_phi = gobjects(nSeta,1);

for k = 1:nSeta
    h_phi(k) = plot( ...
        phi,pphi_all(:,k), ...
        'LineWidth',3, ...
        'Color',color_order(k,:));
end

xlabel('$\varphi$', ...
    'Interpreter','latex', ...
    'FontSize',24, ...
    'FontName','Times New Roman');

ylabel('$p(\varphi)$', ...
    'Interpreter','latex', ...
    'FontSize',24, ...
    'FontName','Times New Roman');

legend( ...
    h_phi,seta_labels, ...
    'Location','northwest', ...
    'FontSize',24, ...
    'FontName','Times New Roman', ...
    'Interpreter','latex', ...
    'Box','off');

ax = gca;
set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',24, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.5, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');

xlim([0,pi]);
ylim([0,1.30*max(pphi_all(:))]);

xticks([0,pi/2,pi]);
xticklabels({'$0$','$\pi/2$','$\pi$'});

grid on;
box on;

save_all_formats(fig_phi,save_dir,'6_guass_Seta_pphi');

%% 7. SAVE
results = struct();
results.parameters = par;
results.settings = cfg;
results.Omega = Omega;
results.xi1 = par.xi1;
results.xi2 = par.xi2;
results.S_eta = S_eta_list;
results.a = a(:);
results.phi = phi(:);
results.da = da;
results.dphi = dphi;
results.pa = pa_all;
results.pphi = pphi_all;
results.joint_pdf = joint_pdf;
results.diagnostics = diagnostics;

save(fullfile(save_dir,cfg.save_mat),'results','-v7.3');

write_marginal_tables_seta( ...
    a,phi,pa_all,pphi_all,S_eta_list, ...
    fullfile(save_dir,cfg.save_pa_csv), ...
    fullfile(save_dir,cfg.save_pphi_csv));

moment_table = table( ...
    S_eta_list(:), ...
    [diagnostics.Ea].', ...
    [diagnostics.Ea2].', ...
    [diagnostics.Var_a].', ...
    [diagnostics.mode_a].', ...
    [diagnostics.mode_phi].', ...
    'VariableNames', ...
    {'S_eta','E_a','E_a2','Var_a','mode_a','mode_phi'});

writetable( ...
    moment_table, ...
    fullfile(save_dir,cfg.save_moments_csv));

fprintf('\n====================================================\n');
fprintf('FILES SAVED IN: %s\n',save_dir);
fprintf('====================================================\n');
fprintf('%s\n',cfg.save_mat);
fprintf('%s\n',cfg.save_pa_csv);
fprintf('%s\n',cfg.save_pphi_csv);
fprintf('%s\n',cfg.save_moments_csv);

if cfg.quick_debug
    fprintf(['\nNOTE: QUICK mode uses weak reset regularization. ', ...
        'Use the curves to compare peak locations and peak topology, ', ...
        'not final stationary branch occupation probabilities.\n']);
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
function write_marginal_tables_seta(a,phi,pa,pphi,S_eta_list,pa_file,pphi_file)

var_names = cell(1,numel(S_eta_list)+1);
var_names{1} = 'a';
for k = 1:numel(S_eta_list)
    var_names{k+1} = matlab.lang.makeValidName( ...
        sprintf('p_a_Seta_%.8f',S_eta_list(k)));
end
Tpa = array2table([a(:),pa],'VariableNames',var_names);
writetable(Tpa,pa_file);

var_names_phi = cell(1,numel(S_eta_list)+1);
var_names_phi{1} = 'phi';
for k = 1:numel(S_eta_list)
    var_names_phi{k+1} = matlab.lang.makeValidName( ...
        sprintf('p_phi_Seta_%.8f',S_eta_list(k)));
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
    'xi1',NaN, ...
    'xi2',NaN, ...
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


function save_all_formats(fig_handle,save_dir,base_name)
% Save JPG, PNG, PDF, and editable MATLAB FIG files.
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

jpg_file = fullfile(save_dir,[base_name,'.jpg']);
png_file = fullfile(save_dir,[base_name,'.png']);
pdf_file = fullfile(save_dir,[base_name,'.pdf']);
fig_file = fullfile(save_dir,[base_name,'.fig']);

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
