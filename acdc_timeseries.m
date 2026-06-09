function [timestamps, J_out_values, converged_flags, varargout] = acdc_timeseries(timestamp_data, concentration_data, temperature_data, varargin)
%ACDC_TIMESERIES_CUSTOM Time-varying ACDC simulation with fixed step time.
%
% Key behavior in this version:
%   1) Fixed simulation time per timestep
%      - step 1 uses FirstStepTime (near steady-state spin-up)
%      - step >=2 uses StepTime
%   2) No adaptive time-doubling loop
%   3) Optional dynamic get_cs.m update from SMPS PSD
%   4) Saves variables.mat every timestep for ACDCsource restart
%   5) Tracks cluster concentration in size bin (default 1.5-1.8 nm)
%
% Inputs:
%   timestamp_data      datetime/duration/numeric (n x 1)
%   concentration_data  [NH3, H2SO4] in cm^-3
%   temperature_data    K
%
% Name-value options:
%   'FirstStepTime'         default 100   (s)
%   'StepTime'              default 1     (s)
%   'SourcesIn'             default 'sources.txt'
%   'SourcesOut'            default 'source_out.txt'
%   'VariablesFile'         default 'variables.mat'
%   'BinRangeNm'            default [1.5 1.8]
%   'UseMobilityDiameter'   default true
%   'DynamicCS'             default false
%   'SMPSData'              struct with fields: times, diam_nm, dndlog
%   'SMPSLogBase'           'ln' (default) or 'log10'
%   'PressurePa'            scalar or n-vector, default 101325
%   'RhoAero'               kg m^-3, default 1500
%   'GetCSPath'             default 'get_cs.m'
%   'Verbose'               default true
%
% Outputs:
%   timestamps          seconds from start
%   J_out_values        cm^-3 s^-1
%   converged_flags     driver convergence flags
%   details (optional)  struct with fields:
%     .C_bin_prev_sum, .C_bin_sum, .C_bin_prev_each, .C_bin_each, .bin_info

p = inputParser;
p.addParameter('FirstStepTime', 10000000000);
p.addParameter('StepTime', 5);
p.addParameter('SourcesIn', 'sources.txt');
p.addParameter('SourcesOut', 'source_out.txt');
p.addParameter('VariablesFile', 'variables.mat');
p.addParameter('BinRangeNm', [1.406 1.742]);
p.addParameter('UseMobilityDiameter', true);
p.addParameter('DynamicCS', false);
p.addParameter('SMPSData', struct());
p.addParameter('SMPSLogBase', 'ln');
p.addParameter('PressurePa', 101325);
p.addParameter('RhoAero', 1500);
p.addParameter('GetCSPath', 'get_cs.m');
p.addParameter('Verbose', true);
p.parse(varargin{:});
opt = p.Results;

%% Time vectors
abs_times = local_to_seconds(timestamp_data);
if isempty(abs_times)
    error('timestamp_data is empty.');
end
timestamps = abs_times - abs_times(1);
timestamps = timestamps(:);

temperature_data = temperature_data(:);
if size(concentration_data,1) ~= numel(timestamps)
    concentration_data = concentration_data.';
end
if size(concentration_data,1) ~= numel(timestamps) || size(concentration_data,2) < 2
    error('concentration_data must be [n x 2] for [NH3, H2SO4].');
end
if numel(temperature_data) ~= numel(timestamps)
    error('temperature_data length must match timestamps.');
end

n_points = numel(timestamps);

if opt.Verbose
    fprintf('=== ACDC Time-Varying Simulation (Fixed Step Time) ===\n');
    fprintf('Points: %d\n', n_points);
    fprintf('Step 1 time: %.3g s | Other steps: %.3g s\n', opt.FirstStepTime, opt.StepTime);
    fprintf('Bin range: [%.3f, %.3f] nm\n', opt.BinRangeNm(1), opt.BinRangeNm(2));
    fprintf('DynamicCS: %d\n\n', logical(opt.DynamicCS));
end

%% Cluster geometry from driver definitions
[clust_names, dp_nm, m_amu] = get_acdc_cluster_props(opt.UseMobilityDiameter);
ncl = numel(dp_nm);

bin_idx = find(dp_nm >= opt.BinRangeNm(1) & dp_nm <= opt.BinRangeNm(2));
if isempty(bin_idx)
    warning('No cluster diameters fall inside [%.3f, %.3f] nm.', opt.BinRangeNm(1), opt.BinRangeNm(2));
end

%% Pressure vector
if isscalar(opt.PressurePa)
    P_vec = repmat(double(opt.PressurePa), n_points, 1);
else
    P_vec = double(opt.PressurePa(:));
    if numel(P_vec) ~= n_points
        error('PressurePa must be scalar or length n_points.');
    end
end

%% Optional SMPS setup
use_dynamic_cs = logical(opt.DynamicCS);
if use_dynamic_cs
    if ~isfield(opt.SMPSData, 'times') || ~isfield(opt.SMPSData, 'diam_nm') || ~isfield(opt.SMPSData, 'dndlog')
        error('DynamicCS=true requires SMPSData struct with fields times, diam_nm, dndlog.');
    end

    smps_t = local_to_seconds(opt.SMPSData.times);
    smps_dp_nm = double(opt.SMPSData.diam_nm(:)).';
    smps_dnd = double(opt.SMPSData.dndlog);

    if size(smps_dnd,1) ~= numel(smps_t) && size(smps_dnd,2) == numel(smps_t)
        smps_dnd = smps_dnd.';
    end
    if size(smps_dnd,1) ~= numel(smps_t)
        error('SMPS dndlog rows must match SMPS times.');
    end
    if size(smps_dnd,2) ~= numel(smps_dp_nm)
        error('SMPS dndlog columns must match number of SMPS diameter bins.');
    end

    % Treat missing values as zero (requested behavior)
    smps_dnd(~isfinite(smps_dnd)) = 0;
    smps_dnd = max(smps_dnd, 0);

    % sort by time for interpolation
    [smps_t, isrt] = sort(smps_t);
    smps_dnd = smps_dnd(isrt, :);

    % interpolate all bins to model times
    dnd_interp = zeros(n_points, size(smps_dnd,2));
    for j = 1:size(smps_dnd,2)
        dnd_interp(:,j) = interp1(smps_t, smps_dnd(:,j), abs_times, 'linear', 'extrap');
    end
    dnd_interp(~isfinite(dnd_interp)) = 0;
    dnd_interp = max(dnd_interp, 0);
else
    dnd_interp = [];
    smps_dp_nm = [];
end

%% Outputs
J_out_values = NaN(n_points,1);
converged_flags = zeros(n_points,1);

C_bin_prev_sum = NaN(n_points,1);
C_bin_sum = NaN(n_points,1);
C_bin_prev_each = NaN(n_points, numel(bin_idx));
C_bin_each = NaN(n_points, numel(bin_idx));

%% Time loop
for i = 1:n_points
    verbose_i = opt.Verbose && (i <= 10 || mod(i,50)==0 || i > n_points-5);

    if verbose_i
        fprintf('Step %d/%d (t=%.3f h): ', i, n_points, timestamps(i)/3600);
    end

    try
        nh3_conc = concentration_data(i,1);
        h2so4_conc = concentration_data(i,2);
        T_K = temperature_data(i);

        % Read previous-step concentrations from variables.mat (before ACDCsource overwrite)
        [prev_sum, prev_each] = extract_prev_bin_from_variables(opt.VariablesFile, bin_idx, ncl);
        C_bin_prev_sum(i) = prev_sum;
        if ~isempty(prev_each)
            C_bin_prev_each(i,1:numel(prev_each)) = prev_each(:).';
        end

        % Update sources/initials from previous run
        ACDCsource(nh3_conc, h2so4_conc);

        % Optional dynamic CS: overwrite get_cs.m for this step
        if use_dynamic_cs
            cs_vec = calc_cs_vec_from_smps_row(smps_dp_nm, dnd_interp(i,:), T_K, ...
                'P_Pa', P_vec(i), ...
                'rho_aero', opt.RhoAero, ...
                'log_base', opt.SMPSLogBase, ...
                'dp_clust_nm', dp_nm, ...
                'm_clust_amu', m_amu);
            write_get_cs_dynamic(cs_vec, opt.GetCSPath, clust_names);
        end

        % Fixed simulation time per step
        if i == 1
            sim_time = opt.FirstStepTime;  % near-steady-state spin-up
        else
            sim_time = opt.StepTime;       % fixed short step afterward
        end

        % Single run (no adaptive loop)
        [C, T, conv_status, clust, Cf, labels_ch, clust_flux, J_out, flux] = ... %#ok<ASGLU>
            driver_acdc(sim_time, 'Sources_in', opt.SourcesIn, 'Sources_out', opt.SourcesOut);

        if ~isempty(J_out)
            J_out_values(i) = J_out(end);
        else
            J_out_values(i) = NaN;
        end
        converged_flags(i) = conv_status;

        % Current-step bin concentration from final C row
        if ~isempty(C)
            c_last = C(end,:);
            n_use = min(numel(c_last), ncl);
            c_last = c_last(1:n_use);

            idx_valid = bin_idx(bin_idx <= n_use);
            if ~isempty(idx_valid)
                c_each = c_last(idx_valid);
                C_bin_each(i,1:numel(c_each)) = c_each;
                C_bin_sum(i) = sum(c_each);
            else
                C_bin_sum(i) = 0;
            end
        end

        % Save restart variables for next timestep (ACDCsource reads this)
        try
            save(opt.VariablesFile, 'C', 'flux', 'clust', 'Cf', 'labels_ch', 'clust_flux', 'J_out', 'T');
        catch
            % If full save fails, at least keep essential fields
            save(opt.VariablesFile, 'C', 'flux', 'clust');
        end

        if verbose_i
            fprintf('sim_time=%.3g s, J=%.3e, conv=%d, Cbin(prev/curr)=%.3e / %.3e cm^-3\n', ...
                sim_time, J_out_values(i), converged_flags(i), C_bin_prev_sum(i), C_bin_sum(i));
        end

    catch ME
        converged_flags(i) = -999;
        if verbose_i
            fprintf('ERROR: %s\n', ME.message);
        end
    end
end

%% Summary
if opt.Verbose
    valid = isfinite(J_out_values);
    fprintf('\n=== Summary ===\n');
    fprintf('Valid J points: %d / %d\n', sum(valid), n_points);
    if any(valid)
        fprintf('J min/max/mean: %.3e / %.3e / %.3e cm^-3 s^-1\n', min(J_out_values(valid)), max(J_out_values(valid)), mean(J_out_values(valid)));
    end
end

%% Optional detail output
details = struct();
details.C_bin_prev_sum = C_bin_prev_sum;
details.C_bin_sum = C_bin_sum;
details.C_bin_prev_each = C_bin_prev_each;
details.C_bin_each = C_bin_each;
details.bin_info = struct( ...
    'indices', bin_idx(:), ...
    'names', {clust_names(bin_idx)}, ...
    'diameters_nm', dp_nm(bin_idx), ...
    'range_nm', opt.BinRangeNm, ...
    'use_mobility', logical(opt.UseMobilityDiameter));

if nargout > 3
    varargout{1} = details;
end

end

%% ========================= helpers =========================
function [sum_bin, each_bin] = extract_prev_bin_from_variables(vars_file, bin_idx, ncl)
sum_bin = NaN;
each_bin = [];

if isempty(bin_idx)
    sum_bin = 0;
    each_bin = [];
    return;
end

if exist(vars_file, 'file') ~= 2
    return;
end

S = load(vars_file);
if ~isfield(S, 'C') || isempty(S.C)
    return;
end

C = S.C;
if size(C,1) < 1
    return;
end

c_last = C(end,:);
n_use = min(numel(c_last), ncl);
c_last = c_last(1:n_use);
idx_valid = bin_idx(bin_idx <= n_use);

if isempty(idx_valid)
    sum_bin = 0;
    each_bin = [];
else
    each_bin = c_last(idx_valid);
    sum_bin = sum(each_bin);
end
end

function [clust, dp_nm, m_amu] = get_acdc_cluster_props(use_mobility)
if nargin < 1, use_mobility = true; end

clust = {'1sa' '1dma' '1dma1sa' '1dma2sa' '2dma1sa' '2dma2sa' '2dma3sa' ...
         '3dma2sa' '3dma3sa' '3dma4sa' '4dma3sa' '4dma4sa' '4dma5sa' ...
         '5dma4sa' '5dma5sa' '5dma6sa' '6dma5sa' '6dma6sa' '6dma7sa' ...
         '7dma6sa' '7dma7sa' '7dma8sa' '8dma7sa' '8dma8sa' '8dma9sa' ...
         '9dma8sa' '9dma9sa' '9dma10sa' '10dma9sa' '10dma10sa' ...
         '10dma11sa' '11dma10sa' '11dma11sa'};

% from driver_acdc.m
diameters = [0.55 0.59 0.72 0.82 0.84 0.91 0.98 0.99 1.04 1.09 ...
             1.11 1.15 1.19 1.20 1.24 1.27 1.28 1.32 1.35 1.36 ...
             1.39 1.41 1.42 1.45 1.48 1.48 1.51 1.53 1.54 1.56 ...
             1.58 1.59 1.61];

mobility_diameters = [0.85 0.89 1.02 1.12 1.14 1.21 1.28 1.29 1.34 1.39 ...
                      1.41 1.45 1.49 1.50 1.54 1.57 1.58 1.62 1.65 1.66 ...
                      1.69 1.71 1.72 1.75 1.78 1.78 1.81 1.83 1.84 1.86 ...
                      1.88 1.89 1.91];

m_amu = [98.08 45.08 143.16 241.24 188.24 286.32 384.40 331.40 429.48 527.56 ...
         474.56 572.64 670.72 617.72 715.80 813.88 760.88 858.96 957.04 904.04 ...
         1002.12 1100.20 1047.20 1145.28 1243.36 1190.36 1288.44 1386.52 1333.52 ...
         1431.60 1529.68 1476.68 1574.76];

if use_mobility
    dp_nm = mobility_diameters(:);
else
    dp_nm = diameters(:);
end
m_amu = m_amu(:);
end

function s = local_to_seconds(t)
if isdatetime(t)
    tt = t(:);
    if isempty(tt.TimeZone)
        tt.TimeZone = 'UTC';
    end
    s = posixtime(tt);
elseif isduration(t)
    s = seconds(t(:));
else
    s = double(t(:));
end
end

function [cs_vec, c_aero_bin] = calc_cs_vec_from_smps_row(diam_nm, dndlog_row, T_K, varargin)
% Follows Fortran get_cs_aero implementation.
p = inputParser;
p.addParameter('P_Pa', 101325);
p.addParameter('rho_aero', 1500);
p.addParameter('log_base', 'ln');
p.addParameter('dp_clust_nm', []);
p.addParameter('m_clust_amu', []);
p.parse(varargin{:});
opt = p.Results;

diam_nm = double(diam_nm(:)).';
dndlog_row = double(dndlog_row(:)).';
if numel(diam_nm) ~= numel(dndlog_row)
    error('diam_nm and dndlog_row must have same length.');
end

if isempty(opt.dp_clust_nm) || isempty(opt.m_clust_amu)
    [~, dp_clust_nm, m_clust_amu] = get_acdc_cluster_props(true);
else
    dp_clust_nm = opt.dp_clust_nm(:);
    m_clust_amu = opt.m_clust_amu(:);
end

ncl = numel(dp_clust_nm);

[dp_aero_m, idx] = sort(diam_nm*1e-9, 'ascend');
dndlog_row = dndlog_row(idx);
% --- ADD THIS FIX TO PREVENT DOUBLE COUNTING ---
% Find the maximum size of the ACDC clusters
%max_clust_m = max(dp_clust_m);

% Keep only SMPS bins that are strictly larger than the ACDC clusters
%valid_idx = dp_aero_m > max_clust_m;
valid_idx = dp_aero_m > 5e-8;
dp_aero_m = dp_aero_m(valid_idx);
dndlog_row = dndlog_row(valid_idx);
% -----------------------------------------------


dndlog_row(~isfinite(dndlog_row)) = 0;
dndlog_row = max(dndlog_row,0);

% dln = local_dln_from_centers(dp_aero_m);
dln = 0.116 ;

switch lower(string(opt.log_base))
    case "ln"
        c_aero_bin = dndlog_row .* dln * 1e6; % cm^-3 -> m^-3
    case "log10"
        c_aero_bin = dndlog_row .* (dln/log(10)) * 1e6;
    otherwise
        error('log_base must be ''ln'' or ''log10''.');
end
c_aero_bin = c_aero_bin(:);

dp_clust_m  = dp_clust_nm * 1e-9;
mp_clust_kg = m_clust_amu * 1e-3 / 6.02214179e23;
mp_aero_kg  = opt.rho_aero * (pi/6) * (dp_aero_m(:).^3);

m_air = 28.97e-3;
Rg = 8.3145;
kB = 1.3806504e-23;
P = opt.P_Pa;
T = T_K;

dp_all = [dp_clust_m; dp_aero_m(:)];
m_all  = [mp_clust_kg; mp_aero_kg];

mu = 1.8e-5*(T/298)^0.85;
lambda = 2*mu/(P*sqrt(8*m_air/(pi*Rg*T)));

Cc = 1 + 2*lambda./dp_all .* (1.257 + 0.4*exp(-1.1*dp_all/(2*lambda)));
Diff = Cc.*kB*T ./ (3*pi*mu*dp_all);

veloc = sqrt(8*kB*T ./ (pi*m_all));
mfp = 8*Diff ./ (pi*veloc);

% Create a radius variable
r_all = dp_all / 2;

% Calculate g using RADIUS, not diameter
g = (1./(3*r_all.*mfp)) .* ((r_all+mfp).^3 - (r_all.^2 + mfp.^2).^(3/2)) - r_all;
g = max(g,0);

na = numel(dp_aero_m);
cs_vec = zeros(ncl,1);

for i = 1:ncl
    acc = 0;
    for j = 1:na
        jj = ncl + j;
        dij = dp_all(i) + dp_all(jj);

        fs_corr = 1 / ( ...
            dij/(dij + 2*sqrt(g(i)^2 + g(jj)^2)) + ...
            8*(Diff(i)+Diff(jj)) / (sqrt(veloc(i)^2 + veloc(jj)^2)*dij) );

        alpha = 0.5; % Example: Only 50% of collisions result in sticking
	beta_ij = alpha*2*pi*dij*(Diff(i)+Diff(jj))*fs_corr;
        acc = acc + beta_ij*c_aero_bin(j);
    end
    cs_vec(i) = acc;
end
end

function dln = local_dln_from_centers(dp_centers_m)
dp = dp_centers_m(:);
n = numel(dp);
if n < 2
    error('Need >= 2 diameter bins.');
end
lnD = log(dp);
lnE = zeros(n+1,1);
lnE(2:n) = 0.5*(lnD(1:n-1)+lnD(2:n));
lnE(1) = lnD(1) - 0.5*(lnD(2)-lnD(1));
lnE(n+1) = lnD(n) + 0.5*(lnD(n)-lnD(n-1));
dln = diff(lnE).';
end

function write_get_cs_dynamic(cs_vec, get_cs_path, clust_names)
if nargin < 2 || isempty(get_cs_path)
    get_cs_path = 'get_cs.m';
end
if nargin < 3 || isempty(clust_names)
    [clust_names, ~, ~] = get_acdc_cluster_props(true);
end

fid = fopen(get_cs_path, 'w');
if fid < 0
    error('Cannot write %s', get_cs_path);
end
cleanObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '%% Auto-generated by acdc_timeseries (dynamic SMPS CS)\n');
for i = 1:numel(cs_vec)
    fprintf(fid, 'CS(%d) = %.16e;\t%% %s\n', i, cs_vec(i), clust_names{i});
end
end
