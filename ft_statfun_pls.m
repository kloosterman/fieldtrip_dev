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
  %			1. Mean-Centering Task PLS
  %			2. Non-Rotated Task PLS
  %			3. Regular Behavior PLS
  %			4. Regular Multiblock PLS
  %			5. Non-Rotated Behavior PLS
  %			6. Non-Rotated Multiblock PLS
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
if ~isfield(cfg, 'zscorescores')
  cfg.zscorescores = 'no';
end

ngroups = length(cfg.num_subj_lst);
ncond   = cfg.num_cond;

datamat_lst = cell(1, ngroups);

% Map columns of 'data' to groups (groups sized as num_subj * ncond)
indexVector = groupSizesToIndices(cfg.num_subj_lst * ncond);

% If doing Behavior PLS, ensure design has rows = (subjects*conditions across all groups),
% and columns = behavior measures
if cfg.pls_method == 3
    if size(design,1) < size(design,2)   % make rows = observations
        designT = design.';
    else
        designT = design;
    end
    behav_blocks = cell(1, ngroups);  % temp, then we'll vertcat into 2-D
end

for i = 1:ngroups
    % --------- brain/data block for group i
    groupdat = data(:, indexVector==i).';
    groupdat_incond = zeros(size(groupdat));

    % within-group condition index (subjects repeated per condition)
    indexVector_cond = groupSizesToIndices(repmat(cfg.num_subj_lst(i),1,ncond));

    for icond = 1:ncond
        conddat = groupdat(indexVector_cond==icond, :);
        groupdat_incond(icond:ncond:end, :) = conddat;
    end
    datamat_lst{i} = groupdat_incond;

    % --------- behavior/design block for group i (matches row order!)
    if cfg.pls_method == 3
        groupdes = designT(indexVector==i, :);         % pick rows for this group
        groupdes_incond = zeros(size(groupdes));
        for icond = 1:ncond
            conddes = groupdes(indexVector_cond==icond, :);
            groupdes_incond(icond:ncond:end, :) = conddes;
        end
        behav_blocks{i} = groupdes_incond;             % rows align with datamat_lst{i}
    end
end

% --------- PLS options
pls_options = [];
pls_options.num_perm     = cfg.num_perm;
pls_options.cormode      = cfg.cormode;
pls_options.num_boot     = cfg.num_boot;
pls_options.method       = cfg.pls_method;
pls_options.num_subj_lst = cfg.num_subj_lst;

% Behavior PLS / Multiblock: single 2-D matrix, rows == sum of rows in datamat_lst
if cfg.pls_method == 3
    pls_options.stacked_behavdata = vertcat(behav_blocks{:});  % [sum_rows x #behav_measures]
    % sanity check (optional):
    % assert(size(pls_options.stacked_behavdata,1) == sum(cellfun(@(x)size(x,1), datamat_lst)), ...
    %    'Row mismatch between behavior and data blocks');
end

% Check for any additional PLS-specific configuration in cfg
fields = fieldnames(cfg);
for i = 1:length(fields)
  if startsWith(fields{i}, 'pls_')
    option_name = strrep(fields{i}, 'pls_', '');
    pls_options.(option_name) = cfg.(fields{i});
  end
end

%% Call pls_analysis with num_cond as input 3
[results] = pls_analysis(datamat_lst, cfg.num_subj_lst, cfg.num_cond, pls_options);
%%

% Parse PLS results into FieldTrip-compatible structure
stat = struct();
stat.latent = results.usc;
% stat.brainscores = results.usc(:,1);
stat.brainscores = cell(1,ngroups);
stat.behavscores = cell(1,ngroups);
if strcmp(cfg.zscorescores, 'yes')
  results.usc = zscore(results.usc);
  results.vsc = zscore(results.vsc);
end
for i=1:ngroups
  stat.brainscores{i} = results.usc(indexVector==i,1);
  stat.brainscores{i} = reshape(stat.brainscores{i}, [], cfg.num_cond);
  stat.behavscores{i} = results.vsc(indexVector==i,1);
  stat.behavscores{i} = reshape(stat.behavscores{i}, [], cfg.num_cond);
end
stat.results = results;
stat.cfg = cfg;

if isfield(results, 'boot_result')
  stat.boot_res = results.boot_result;
  stat.stat = results.boot_result.compare_u(:,1); % bootstrap ratios
end

if isfield(results, 'perm_result')
  for i = 1:length(results.perm_result.sprob)
    stat.posclusters(i).prob = results.perm_result.sprob(i);
  end
  stat.perm_res = results.perm_result;
end

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
