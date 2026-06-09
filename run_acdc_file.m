function run_acdc_file(filename, varargin)
%RUN_ACDC_FILE Run ACDC time-series from main forcing CSV.
%
% Drop-in replacement with:
%   - fixed step-time simulation control (step1 vs later steps)
%   - optional dynamic CS from SMPS file
%   - SMPS missing/empty values treated as zero
%
% Usage:
%   run_acdc_file('processed.csv')
%   run_acdc_file('processed.csv', 'FirstStepTime', 100, 'StepTime', 1)
%   run_acdc_file('processed.csv', 'DynamicCS', true, 'SMPSFile', 'smps.csv')

p = inputParser;
p.addRequired('filename', @(x)ischar(x) || isstring(x));
p.addParameter('FirstStepTime', 1000000000);
p.addParameter('StepTime', 5);
p.addParameter('BinRangeNm', [1.5 1.8]);
p.addParameter('UseMobilityDiameter', true);
p.addParameter('AmmoniaConc', 7.0e7); % cm^-3 constant NH3

% Dynamic-CS options
p.addParameter('DynamicCS', false);
p.addParameter('SMPSFile', '');
p.addParameter('SMPSLogBase', 'ln');     % 'ln' or 'log10'
p.addParameter('PressurePa', 101325);
p.addParameter('RhoAero', 1500);
p.addParameter('GetCSPath', 'get_cs.m');

p.parse(filename, varargin{:});
opt = p.Results;
filename = char(opt.filename);

if exist(filename, 'file') ~= 2
    error('File does not exist: %s', filename);
end

fprintf('Processing forcing file: %s\n', filename);

%% Read primary forcing CSV
[timestamps_dt, temperatures, sa_concentrations] = read_main_forcing_csv(filename);

% Constant NH3 (edit here if you have NH3 timeseries)
ammonia_concentrations = opt.AmmoniaConc * ones(size(sa_concentrations));
concentrations = [ammonia_concentrations, sa_concentrations];

fprintf('\nForcing summary:\n');
fprintf('  Points: %d\n', numel(temperatures));
fprintf('  Time span: %.2f h\n', hours(timestamps_dt(end)-timestamps_dt(1)));
fprintf('  Temperature: [%.2f, %.2f] K\n', min(temperatures), max(temperatures));
fprintf('  SA: [%.3e, %.3e] cm^-3\n', min(sa_concentrations), max(sa_concentrations));
fprintf('  NH3 (const): %.3e cm^-3\n', opt.AmmoniaConc);

%% Optional SMPS data (for dynamic CS)
smpsData = struct();
if opt.DynamicCS
    if isempty(opt.SMPSFile)
        error('DynamicCS=true requires ''SMPSFile''.');
    end
    if exist(opt.SMPSFile, 'file') ~= 2
        error('SMPS file not found: %s', opt.SMPSFile);
    end

    [smps_times, smps_diam_nm, smps_dnd] = process_smps_data(opt.SMPSFile);
    smpsData.times = smps_times;
    smpsData.diam_nm = smps_diam_nm;
    smpsData.dndlog = smps_dnd;

    fprintf('\nSMPS summary:\n');
    fprintf('  File: %s\n', opt.SMPSFile);
    fprintf('  Rows: %d | Bins: %d\n', size(smps_dnd,1), size(smps_dnd,2));
    fprintf('  Diameter range: [%.3f, %.3f] nm\n', min(smps_diam_nm), max(smps_diam_nm));
    fprintf('  Missing/empty entries treated as 0\n');
end

%% Run simulation
fprintf('\nRunning ACDC simulation...\n');
[times_s, J_out, converged, details] = acdc_timeseries( ...
    timestamps_dt, concentrations, temperatures, ...
    'FirstStepTime', opt.FirstStepTime, ...
    'StepTime', opt.StepTime, ...
    'BinRangeNm', opt.BinRangeNm, ...
    'UseMobilityDiameter', opt.UseMobilityDiameter, ...
    'DynamicCS', opt.DynamicCS, ...
    'SMPSData', smpsData, ...
    'SMPSLogBase', opt.SMPSLogBase, ...
    'PressurePa', opt.PressurePa, ...
    'RhoAero', opt.RhoAero, ...
    'GetCSPath', opt.GetCSPath);

%% Save outputs
[folder, base, ~] = fileparts(filename);
if isempty(folder), folder = pwd; end

out_csv = fullfile(folder, sprintf('%s_acdc_results.csv', base));
out_mat = fullfile(folder, sprintf('%s_acdc_results.mat', base));

results_table = table( ...
    timestamps_dt, times_s, temperatures, sa_concentrations, ammonia_concentrations, ...
    J_out, converged, details.C_bin_prev_sum, details.C_bin_sum, ...
    'VariableNames', {'DateTime','Time_s','Temperature_K','SA_conc','NH3_conc', ...
                      'J_out','Converged','C_bin_prev_sum','C_bin_sum'});

writetable(results_table, out_csv);
save(out_mat, 'results_table', 'details', 'opt');

fprintf('\nSaved:\n  %s\n  %s\n', out_csv, out_mat);

%% Brief report
valid = isfinite(J_out);
fprintf('\nResult summary:\n');
fprintf('  Valid J points: %d/%d\n', sum(valid), numel(valid));
if any(valid)
    fprintf('  J min/max/mean: %.3e / %.3e / %.3e cm^-3 s^-1\n', min(J_out(valid)), max(J_out(valid)), mean(J_out(valid)));
end

fprintf('Done.\n');
end

%% ========================= helpers =========================
function [timestamps_dt, temperatures, sa_concentrations] = read_main_forcing_csv(filename)
% Flexible CSV read: expects at least timestamp, temperature, SA columns.

try
    T = readtable(filename, 'Delimiter', ',', 'ReadVariableNames', true, ...
        'VariableNamingRule', 'preserve', 'TextType', 'string');
catch ME
    error('Failed to read forcing CSV %s: %s', filename, ME.message);
end

if height(T) == 0
    error('Forcing CSV has no rows: %s', filename);
end

col_names = T.Properties.VariableNames;

% Detect columns by name; fallback to positions 1/2/3
idx_time = find(contains(lower(col_names), 'time') | contains(lower(col_names), 'date'), 1);
idx_temp = find(contains(lower(col_names), 'temp'), 1);
idx_sa   = find(contains(lower(col_names), 'sa') | contains(lower(col_names), 'sulfur'), 1);

if isempty(idx_time), idx_time = 1; end
if isempty(idx_temp), idx_temp = 2; end
if isempty(idx_sa),   idx_sa = 3;   end

raw_time = T{:, idx_time};
temperatures = to_numeric_col(T{:, idx_temp});
sa_concentrations = to_numeric_col(T{:, idx_sa});

% Clean NaNs
tmp_med = median(temperatures(isfinite(temperatures)));
if isempty(tmp_med) || ~isfinite(tmp_med), tmp_med = 298.15; end
temperatures(~isfinite(temperatures)) = tmp_med;
sa_concentrations(~isfinite(sa_concentrations)) = 0;

% Parse datetime
if isdatetime(raw_time)
    timestamps_dt = raw_time;
elseif isduration(raw_time)
    t0 = datetime('now');
    timestamps_dt = t0 + raw_time;
else
    try
        timestamps_dt = datetime(raw_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ssXXX', 'TimeZone', 'UTC');
        timestamps_dt.TimeZone = '';
    catch
        try
            timestamps_dt = datetime(raw_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        catch
            timestamps_dt = datetime(raw_time);
        end
    end
end

timestamps_dt = timestamps_dt(:);
end

function [times, diameters_nm, dndlndp_values] = process_smps_data(file_path)
%PROCESS_SMPS_DATA
% MATLAB equivalent of the user's Python function with robust missing handling.
%
% Expected layout:
%   col 1      : scan start time
%   cols 2..51 : diameter columns (50 bins)
%   cols 52..101: dN/dlnDp columns (50 bins)
%
% Missing/empty values are treated as zero in dN/dlnDp.

T = readtable(file_path, 'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'string');
if width(T) < 3
    error('SMPS file appears invalid (too few columns): %s', file_path);
end

% time column
raw_time = T{:,1};
if isdatetime(raw_time)
    times = raw_time;
else
    try
        times = datetime(raw_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ssXXX', 'TimeZone', 'UTC');
        times.TimeZone = '';
    catch
        try
            times = datetime(raw_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        catch
            times = datetime(raw_time);
        end
    end
end

ncol = width(T);

% match the Python slicing when possible
if ncol >= 101
    diam_idx = 2:51;
    dnd_idx = 52:101;
else
    % fallback: split remaining columns in half
    rem = ncol - 1;
    nbin = floor(rem/2);
    diam_idx = 2:(1+nbin);
    dnd_idx = (2+nbin):(1+2*nbin);
end

diam_row = T{1, diam_idx};
diameters_nm = to_numeric_row(diam_row);

% If diameter entries are empty, infer from valid neighbors is not reliable;
% keep only positive finite bins.
valid_dp = isfinite(diameters_nm) & diameters_nm > 0;
if ~any(valid_dp)
    error('No valid SMPS diameters found in %s', file_path);
end

dnd_raw = T{:, dnd_idx};
dnd = to_numeric_matrix(dnd_raw);

% Requested behavior: empty columns/values => zero
% readtable often maps empty strings to <missing> / NaN; force to zero.
dnd(~isfinite(dnd)) = 0;

% Keep only bins with valid diameters
diameters_nm = diameters_nm(valid_dp);
dndlndp_values = dnd(:, valid_dp);

% Any negative values are nonphysical for number density
% (keep zero floor)
dndlndp_values = max(dndlndp_values, 0);
end

function x = to_numeric_col(v)
if isnumeric(v)
    x = double(v(:));
elseif isstring(v) || iscellstr(v) || iscell(v)
    x = str2double(string(v(:)));
else
    x = double(v(:));
end
end

function x = to_numeric_row(v)
if isnumeric(v)
    x = double(v);
elseif isstring(v) || iscellstr(v) || iscell(v)
    x = str2double(string(v));
else
    x = double(v);
end
x = x(:).';
end

function M = to_numeric_matrix(v)
if isnumeric(v)
    M = double(v);
elseif isstring(v)
    M = str2double(v);
elseif iscell(v)
    M = str2double(string(v));
else
    M = double(v);
end
end
