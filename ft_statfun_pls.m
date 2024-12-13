function stat = ft_statfun_pls(cfg, data, design)
% FT_STATFUN_PLS Wrapper function for Partial Least Squares (PLS) analysis using pls_analysis.m
% 
% This function integrates PLS analysis into FieldTrip's statistical pipeline.
% 
% INPUT:
%   cfg    - Configuration structure with relevant PLS parameters.
%   data   - Data structure with FieldTrip data format.
%   design - Design matrix for the analysis.
% 
% OUTPUT:
%   stat   - Struct with PLS results.
% 
% REQUIREMENTS:
%   pls_analysis.m from Rotman Baycrest software.
% 
% Example usage:
%   cfg.pls_method = 3; % 1 for mean-centering, 2 for non-rotated, 3 for behavioral PLS
%   cfg.statistic = 'ft_statfun_pls';
%   cfg.num_perm = 500;
%   cfg.num_boot = 1000;
% 
% PLS Method Legend:
%   1 - Mean-Centering PLS: Standard PLS method with mean centering of the data.
%   2 - Non-Rotated PLS: An alternative PLS method without rotation.
%   3 - Behavioral PLS: PLS designed for analyzing behavioral data relationships.
% 
% Check if pls_analysis.m is available
if ~exist('pls_analysis', 'file')
    error('pls_analysis.m not found in the MATLAB path. Please add it to your path.');
end

% Parse input data and configuration
if ~isfield(cfg, 'num_perm')
    cfg.num_perm = 500; % Default number of permutations
end
if ~isfield(cfg, 'num_boot')
    cfg.num_boot = 1000; % Default number of bootstrap samples
end
if ~isfield(cfg, 'pls_method')
    cfg.pls_method = 1; % Default to mean-centering PLS
end

% Check for required data.dimord
if ~isfield(data, 'dimord')
    error('The input data must contain a dimord field indicating the data dimension order.');
end

data_dimord = data.dimord;
data_dims = strsplit(data_dimord, '_');

% Assume the data is in data.powspctrm
if ~isfield(data, 'powspctrm')
    error('The data structure must contain a ''powspctrm'' field.');
end

data_field = data.powspctrm;

% Determine the order of dimensions
switch data_dimord
    case 'chan_time_subject'
        data_reshaped = reshape(data_field, size(data_field, 1) * size(data_field, 2), size(data_field, 3))';
    case 'chan_freq_subject'
        data_reshaped = reshape(data_field, size(data_field, 1) * size(data_field, 2), size(data_field, 3))';
    case 'chan_time_freq_subject'
        data_reshaped = reshape(data_field, size(data_field, 1) * size(data_field, 2) * size(data_field, 3), size(data_field, 4))';
    case {'subj_chan_freq_time', 'rpt_chan_freq_time'}
        data_reshaped = reshape(data_field, size(data_field, 2) * size(data_field, 3) * size(data_field, 4), size(data_field, 1));
    otherwise
        error('Unsupported dimord: %s', data_dimord);
end

% Configure PLS options
pls_options.num_perm = cfg.num_perm;
pls_options.num_boot = cfg.num_boot;
pls_options.method = cfg.pls_method; % Directly assign method as a number, as required by pls_analysis

% Check for any additional PLS-specific configuration in cfg
fields = fieldnames(cfg);
for i = 1:length(fields)
    if startsWith(fields{i}, 'pls_')
        option_name = strrep(fields{i}, 'pls_', '');
        pls_options.(option_name) = cfg.(fields{i});
    end
end

[results] = pls_analysis(data_reshaped, design, pls_options);

% Parse PLS results into FieldTrip-compatible structure
stat = struct();
stat.scores = results.scores;
stat.latent = results.latent;
stat.perm_res = results.perm_res;
stat.boot_res = results.boot_res;
stat.cfg = cfg;

end
