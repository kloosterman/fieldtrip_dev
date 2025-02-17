function stat = ft_statfun_pls(cfg, data, design)
% FT_STATFUN_PLS Wrapper function for Partial Least Squares (PLS) analysis using pls_analysis.m
% 
% This function integrates PLS analysis into FieldTrip's statistical pipeline.
% 
% INPUT:
%   cfg    - Configuration structure with relevant PLS parameters.
%   data   - Matrix of data with dimensions observations x subjects.
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
%   cfg.num_cond = 2; % Number of conditions
%   cfg.num_subj_lst = [10, 10]; % Number of subjects per condition
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
if ~isfield(cfg, 'num_cond')
    error('cfg.num_cond must be specified to indicate the number of conditions.');
end
if ~isfield(cfg, 'num_subj_lst')
    error('cfg.num_subj_lst must be specified to indicate the number of subjects per condition.');
end
if ~isfield(cfg, 'cormode')
    cfg.cormode = 0; % Pearson
end
if ~isfield(cfg, 'interaction')
    cfg.interaction = 'no'; 
end

% make cell array of groups
ngroups = length(cfg.num_subj_lst);
datamat_lst = cell(1,ngroups);
indexVector = groupSizesToIndices(cfg.num_subj_lst);
for i=1:ngroups
  datamat_lst{i} =  transpose(data(:, indexVector==i));
  if strcmp(cfg.interaction, 'yes') % augment the datamat with contrast data
    datamat_lst{i} = [datamat_lst{i} cfg.contrast(i) .* datamat_lst{i}]; 
  end
end

% Configure PLS options
pls_options.num_perm = cfg.num_perm;
pls_options.cormode = cfg.cormode;
pls_options.num_boot = cfg.num_boot;
pls_options.method = cfg.pls_method; % Directly assign method as a number, as required by pls_analysis
pls_options.num_subj_lst = cfg.num_subj_lst; % Number of subjects per condition
pls_options.stacked_behavdata = transpose(design); % Add design matrix to PLS options

% Check for any additional PLS-specific configuration in cfg
fields = fieldnames(cfg);
for i = 1:length(fields)
    if startsWith(fields{i}, 'pls_')
        option_name = strrep(fields{i}, 'pls_', '');
        pls_options.(option_name) = cfg.(fields{i});
    end
end

% Call pls_analysis with num_cond as input 3
[results] = pls_analysis(datamat_lst, cfg.num_subj_lst, cfg.num_cond, pls_options);

% Parse PLS results into FieldTrip-compatible structure
stat = struct();
if strcmp(cfg.interaction, 'yes') % augment the datamat with contrast data
  ndatapoints = size(results.boot_result.compare_u,1);
  stat.stat = results.boot_result.compare_u(1:(ndatapoints/2), 1);
  stat.prob = results.boot_result.compare_u((ndatapoints/2)+1:end,1);
else
  stat.stat = results.boot_result.compare_u(:,1);
end

for i = 1:length(results.perm_result.sprob)
  stat.posclusters(i).prob = results.perm_result.sprob(i);
end
stat.latent = results.usc;
stat.perm_res = results.perm_result;
stat.boot_res = results.boot_result;
stat.brainscores = results.usc(:,1);
stat.behavscores = results.vsc(:,1);

stat.results = results;
stat.cfg = cfg;

end


function indexVector = groupSizesToIndices(groupSizes)
    % This function transforms a vector of group sizes into an index vector.
    % Each unique group index is repeated according to its size.

    % Initialize the index vector
    indexVector = [];

    % Loop through each group size and append its indices
    for i = 1:length(groupSizes)
        indexVector = [indexVector, i * ones(1, groupSizes(i))];
    end
end
