addpath('/Users/kloosterman/Documents/GitHub/entropyAge')
addpath('/Users/kloosterman/Documents/GitHub/plscmd')
addpath('/Users/kloosterman/Documents/GitHub/fieldtrip_dev')
addpath('/Users/kloosterman/Documents/GitHub/fieldtrip')

datapath = '/Users/kloosterman/projectdata/EntropyAging';
cd(datapath)

behav = readtable('FoPra_Behavioral_Measures.xlsx');
behavNames = behav.Properties.VariableNames;
behav = behav(behav.Age_Group==1,:);
behav = table2array(behav(:,7:19));
behavNames = behavNames(7:19);
load('young_mse_all.mat', 'young_mse_all')

cfg = [];
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 100;                         % Number of permutation
cfg.num_boot = 100;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 8;                            % 0 is Pearson corr, 8 is Spearman
cfg.design = behav;
cfg.num_cond = 1;                           % Number of conditions
cfg.num_subj_lst = size(young_mse_all.powspctrm,1); % Number of subjects per condition

% Step 3: Compute statistics
stat_mse = ft_freqstatistics(cfg, young_mse_all);

%%
cfg =[]; 
cfg.layout=lay;
cfg.parameter = 'stat';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
stat_mse.mask = double(stat_mse.stat > 3 | stat_mse.stat < -3);
cfg.colormap = cmap;
cfg.maskparameter = 'mask';
ft_multiplotTFR(cfg, stat_mse)

f = figure;
tiledlayout(2,2);
cfg.clussign = 'pos';
cfg.clus2plot = 1;   cfg.integratetype = 'trapz'; % mean or trapz
stat_mse.posclusterslabelmat = stat_mse.mask;
ft_clusterplot3D(cfg, stat_mse)

nexttile; s = scatter(stat_mse.brainscores, stat_mse.behavscores, 'MarkerEdgeColor',[1 1 1],  'MarkerFaceColor', [0 0 0], 'LineWidth',1.0, 'SizeData', 40); axis padded; lsline; box on;
xlabel('EEG entropy'); ylabel('Cognitive performance')
title(sprintf('r = %1.2f', corr(stat_mse.brainscores, stat_mse.behavscores)))
% TODO bar plot corrs for each behav var
nexttile; b=bar(corr(stat_mse.brainscores, -behav)); xticklabels(behavNames);
ylabel('Correlation')

saveas(f, 'behavPLS_MSE', 'pdf');
saveas(f, 'behavPLS_MSE', 'png');




% TODO 2 group PLS

