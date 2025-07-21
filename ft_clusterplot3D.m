function plotsuccess = ft_clusterplot3D(cfg, stat)

% FT_EXAMPLEFUNCTION demonstrates to new developers how a FieldTrip function should look like
%
% NK plot integrated TFR and topo for each cluster
%
% Use as
%
% cfg=[];
% cfg.clus2plot = 1;
% cfg.clussign = 'neg';
% cfg.integratetype = 'mean'; % mean or trapz
% cfg.subplotsize = [2 4]; % 2 rows, 2 topo/TFR couples
% cfg.subplotind = [1 2];
% for itrig = 1:2
%   for ifreq = 1:2
%     ft_clusterplot3D(cfg, stat)
%   end
% end
%
% %   cfg.subplotsize               = layout of subplots ([h w], default [3 5])
%
%   outdata = ft_examplefunction(cfg, indata)
% where indata is <<describe the type of data or where it comes from>>
% and cfg is a configuration structure that should contain
%
% <<note that the cfg list should be indented with two spaces
%
%  cfg.option1    = value, explain the value here (default = something)
%  cfg.option2    = value, describe the value here and if needed
%                   continue here to allow automatic parsing of the help
%
% The configuration can optionally contain
%   cfg.option3   = value, explain it here (default is automatic)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also <<give a list of function names, all in capitals>>

% Here come the Copyrights
%
% Here comes the Revision tag, which is auto-updated by the version control system
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will reset ft_warning and show the function help if nargin==0 and return an error
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar    stat % this reads the input data in case the user specified the cfg.inputfile option
ft_preamble provenance stat % this records the time and memory usage at the beginning of the function
% ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% TODO check stat input etc
% get the options
parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
maskparameter = ft_getopt(cfg, 'mask', 'mask');
clus2plot    = ft_getopt(cfg, 'clus2plot');
clussign    = ft_getopt(cfg, 'clussign');
integratetype    = ft_getopt(cfg, 'integratetype', 'mean');
figposition = ft_getopt(cfg, 'figposition', [100 150 600 100*4]);
subplotsize    = ft_getopt(cfg, 'subplotsize', [2 4]);
subplotind    = ft_getopt(cfg, 'subplotind', [1 2]);
cmap = ft_getopt(cfg, 'colormap', parula);
zlim = ft_getopt(cfg, 'zlim', 'maxabs');
zlimTFR = ft_getopt(cfg, 'zlimTFR', 'maxabs'); % TODO make backward compatible
zlimtopo = ft_getopt(cfg, 'zlimtopo', 'maxabs');
layout = ft_getopt(cfg, 'layout');
titleTFR = ft_getopt(cfg, 'titleTFR', ' ');
mask = ft_getopt(stat, 'mask'); % mask from stat struct!
alpha = ft_getopt(cfg, 'alpha', 1);    % number, highest cluster p-value to be plotted max 1 (default = 0.05)
clrbar= ft_getopt(cfg, 'colorbar', 'no'); 
% prepare the layout, this only has to be done once
tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'scalepos', 'elec', 'grad', 'opto', 'showcallinfo'});
layout = ft_prepare_layout(tmpcfg, stat);
ylabl = ft_getopt(cfg, 'ylabel', 'Frequency/Time scale'); 

if ~isfield(stat, [clussign 'clusters'])
    fprintf('No labelmat present')
    return
end

% return if cluster not significant
clus_pval = stat.([clussign 'clusters'])(clus2plot).prob;
fprintf('Cluster p = %g\n', clus_pval)
if clus_pval > alpha
  fprintf('Not plotting')
  plotsuccess = false;
  return
end

% use either cfg specified mask, or make mask from negclusterslabelmat
% this overrides mask cfg input!
if ~isempty(clus2plot) &&  ~isempty(clussign)
  fprintf('Creating mask from cluster %s%d\n', clussign, clus2plot)
  mask = stat.([clussign 'clusterslabelmat']) == clus2plot;
end
if all(mask(:) == false) % in case mask is empty for some reason
  fprintf('mask is empty, nothing to plot\n')
  plotsuccess = 1;
  return
end
  
  
%% prepare data for TFR
TFRdat = stat.(parameter);

freqTFR = [];
if strcmp(integratetype, 'trapz')
  TFRdat(~mask) = 0; % mask out what's outside of cluster
  freqTFR.powspctrm = trapz(TFRdat,1);
elseif strcmp(integratetype, 'mean')
  TFRdat(~mask) = NaN;
  freqTFR.powspctrm = nanmean(TFRdat,1);
%   TFRdat(isnan(TFRdat))=0;
end
freqTFR.mask = double(any(mask,1));
if all(freqTFR.mask(:) == 1) % maybe a bug that nothing is plotted when allmask==1
  freqTFR.mask = any(mask,1);
end
freqTFR.time = stat.time;
freqTFR.freq = stat.freq;
freqTFR.dimord = 'chan_freq_time';
freqTFR.label = {'avg'};

cfgfreq = [];
cfgfreq.maskparameter = 'mask';
cfgfreq.zlim = zlimTFR;
cfgfreq.colorbar = 'no';
cfgfreq.interactive = 'no';
cfgfreq.title = ' ';
cfgfreq.interpolatenan     = 'no';
cfgfreq.figure = 'gca';

% subplot(subplotsize(1),subplotsize(2), subplotind(1));
% tiledlayout(subplotsize); % set when calling?
nexttile
colormap(gca, cmap); hold on

ft_singleplotTFR(cfgfreq, freqTFR);
% tryout for plotting lines
% freqTFR.dimord = 'chan_freq';
% freqTFR = rmfield(freqTFR, 'time')
% ft_singleplotER(cfgfreq, freqTFR);

box on
ax=gca; hold on
plot([0,0], ax.YLim,'k',[0,0], ax.YLim,'k', 'Linewidth', 0.5);
% ax.Position(3) = length(freqTFR.time)*0.015;
title(sprintf('%s', titleTFR), 'Fontsize', 12)
zlim = get(gca, 'CLim'); % use same CLim for both topo and TFR

c = colorbar;
% c.Position = [0.9, 0.55, 0.02, 0.35];
% 
% c.Position(1) = c.Position(1) * 1;
c.Position(1) = c.Position(1)+0.04;
c.Position(2) = c.Position(2)+0.1;
c.Position(3) = c.Position(3) * 0.5;
c.Position(4) = c.Position(4) * 0.3;
% c.Position(3) = 0.015; % 0.005
% c.Position(4) = 0.35;
c.Box = 'off';
% % if strcmp(clrbar, 'no')
% %     c.Visible = 'off';
% % end
xlabel('Time (s)')
ylabel(ylabl)

%% prepare data topo
chansel = any(mask(:,:),2);

freqtopo = [];
topodat = stat.(parameter);
if strcmp(integratetype, 'trapz')
  topodat(~mask) = 0;
  freqtopo.powspctrm = trapz(topodat(:,:),2);
elseif strcmp(integratetype, 'mean')
  topodat(~mask) = NaN;
  freqtopo.powspctrm = squeeze(nanmean(topodat(:,:),2));
end
freqtopo.label = stat.label; % chlabel(LR_subtract_mat(:,2));
freqtopo.dimord = 'chan';
% freqtopo.mask = chansel; % mask cluster, but gives sharp edges
% freqtopo.mask = true(size(chansel)); % no mask

cfgtopo = [];
cfgtopo.layout = layout;
cfgtopo.showcallinfo = 'no';
cfgtopo.feedback = 'no';
% % add neighbours of cluster chans to the cluster, to mask without sharp edges
% % todo spatial upsampling?
% cfg0_neighb = [];
% cfg0_neighb.elec  = 'standard_1020.elc';
% % cfg0_neighb.method    = 'template';
% % cfg0_neighb.template  = 'elec1010_neighb.mat';
% cfg0_neighb.method    = 'distance';
% cfg0_neighb.neighbourdist = 20;
% cfg0_neighb.feedback = 'no';
% neighbours       = ft_prepare_neighbours(cfg0_neighb);
% chans_incnb = unique(vertcat(neighbours.neighblabel)); % get neighbours of chans in cluster
% freqtopo.mask = ismember(freqtopo.label, [freqtopo.label(chansel); chans_incnb]);
% % end TODO

% plot topo
cfgtopo.comment = 'no';
cfgtopo.marker = 'off';
cfgtopo.shading = 'flat';
cfgtopo.style = 'straight_imsat'; %both  straight
cfgtopo.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
cfgtopo.markersize = 3;
cfgtopo.highlight = 'off';
cfgtopo.highlightsymbol = '.';
cfgtopo.highlightsize = 20;
cfgtopo.zlim = zlimtopo;
cfgtopo.highlightchannel = stat.label(chansel);
cfgtopo.parameter = 'powspctrm'; % also when plotting rho
cfgtopo.maskparameter = 'mask';
cfgtopo.interactive = 'no';
cfgtopo.colormap = cmap; 
% cfgtopo.interpolatenan = 'no';
cfgtopo.interpolation = 'v4';
cfgtopo.figure = 'gca';

% subplot(subplotsize(1),subplotsize(2), subplotind(2));
ax = nexttile;
hold on
ft_topoplotTFR(cfgtopo, freqtopo);
% colormap(gca, cmap); 

%plot sensors with size determined by weight
hold on
weights = double(mask);
weights = sum(weights(:,:),2);
weights = weights / max(weights);
%         weights = weights / mean(weights(weights > 0)); % weights centered around 1
%         weights = weights / sum(weights); % make weights sum to 1

[selchan, sellay] = match_str(stat.label, cfgtopo.layout.label);
chanX      = cfgtopo.layout.pos(sellay, 1);
chanY      = cfgtopo.layout.pos(sellay, 2);
chanWidth  = cfgtopo.layout.width(sellay);
chanHeight = cfgtopo.layout.height(sellay);
chanLabel  = cfgtopo.layout.label(sellay);

for isens = find(weights > 0)'
    plot(chanX(isens,1), chanY(isens,1), 'marker', 'o', 'color', 'k', 'markersize', 4 * weights(isens), 'linestyle', 'none', 'LineWidth', 0.75)
%   plot(cfgtopo.layout.pos(isens,1), cfgtopo.layout.pos(isens,2), 'marker', 'o', 'color', 'k', 'markersize', 4 * weights(isens), 'linestyle', 'none', 'LineWidth', 0.75)
  %                                 plot(senspos(isens,1), senspos(isens,2), 'marker', 'x', 'color', 'k', 'markersize', 10, 'linestyle', 'none')
end
% title(sprintf('%s cluster %d\np = %1.3f', clussign, clus2plot, stat.([clussign 'clusters'])(clus2plot).prob), 'FontWeight', 'normal')
title(sprintf('%s cluster %d\np = %1.3f', clussign, clus2plot, clus_pval ), 'FontWeight', 'normal')

plotsuccess = true; % output argument

c = colorbar;
% c.Position(1) = c.Position(1) * 1.1;
% c.Position(3) = c.Position(3) * 0.5;
% c.Position(4) = c.Position(4) * 0.3;

% c.Position(1) = c.Position(1)+0.075;
% c.Position(2) = c.Position(2)+0.1;
% c.Position(3) = 0.015;
% c.Position(4) = 0.35;
c.Position(1) = c.Position(1)+0.04;
c.Position(2) = c.Position(2)+0.1;
c.Position(3) = c.Position(3) * 0.5;
c.Position(4) = c.Position(4) * 0.3;

c.Box = 'off';
% if strcmp(clrbar, 'no')
%     c.Visible = 'off';
% end

% if length(c.Ticks) == 5
%   c.Ticks = c.Ticks([1 3 5]);
% end
% if plorder(iplot) == 1
%   t = text(-0.6 , 0.7 , panelind(irow), 'Fontsize', 12, 'FontWeight', 'bold');
% end


% the ft_postamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

% ft_postamble debug               % this clears the onCleanup function used for debugging in case of an error
% ft_postamble trackconfig         % this converts the config object back into a struct and can report on the unused fields
% ft_postamble previous   stat   % this copies the stat.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
% ft_postamble provenance dataout  % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
% ft_postamble history    dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
% ft_postamble savevar    dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
