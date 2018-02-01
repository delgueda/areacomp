%% CLK part
% taken from plot_click_contextcomp.m

cd ~/code/mycode/psth
dbopen;
baphy_set_path;

disp('loading cell_list... this might take several minutes');

load ('Data/CLK-PFC.mat'); cell_list(:,1)=[];
disp('Done!');

do_clt = 0;
alfa = 0.05;
ncells = length(cell_list);
passives=0;
actives =0;

for i = 1 : ncells
  passives = passives + sum(cell_list{3,i} == 'p');
  actives  = actives  + sum(cell_list{3,i} == 'a');
end

% durations to plot
respdur = 3.6;
presil  = 0.2;

% psth parameters
% psth_sr   = 20;
% smoothwin = 30; % sliding window
options.rasterfs = 30;
options.psth = 0;
options.includeprestim = 1;

psth_sr = options.rasterfs;

% other definitions
letter  = ['abcdefghijklm'];
minisol = 79.5;   % use only cells with isolation >= 79.5%

for i = 1 : ncells
  for j = 1 : length(cell_list{2, i})
    
    cell_list{2,i}{j}.SortInfo.Isolation = [];
    tmpquery = mysql(['SELECT isolation FROM gSingleRaw WHERE cellid="' cell_list{1,i} '" AND rawid="' num2str(cell_list{2,i}{j}.ID) '"']);
    if isempty(tmpquery) || isempty(tmpquery(1).isolation)
      tmpquery(1).isolation = 0; % assign zero isolation to those that have no data (discard)
    end
    
    cell_list{2,i}{j}.SortInfo.Isolation = tmpquery(1).isolation;
    
  end
end

rclk1=zeros(actives, (respdur*psth_sr)+1) .* NaN; %+2 for dPEG A1, +1 for VPr
rclk2=zeros(actives, (respdur*psth_sr)+1) .* NaN;
rclt =zeros(actives, (respdur*psth_sr)+1) .* NaN;

tclk1=zeros(actives, (respdur*psth_sr)+1) .* NaN;
tclk2=zeros(actives, (respdur*psth_sr)+1) .* NaN;
tclt =zeros(actives, (respdur*psth_sr)+1) .* NaN;

cellid_used = cell(actives,1);  % to save the cellids selected and saved
parameters  = cell(actives,1);  % to save exptparams of the active file

% main loop
% will assume that always before an active there was a passive
cnt = 1;
sel_cell_num = 0;
clt_cnt = 0;

for i = 1 : ncells
  actidx = find(cell_list{3,i} == 'a');
  if actidx(1) ~= 1 
    preidx = actidx - 1;
  else
    preidx = actidx + 1;
  end
  
  if length(cell_list{3,i}) == 1
    preidx = [];
    
  end
  


  % prepare some options
  if strcmp(cell_list{1,i}(end-1), '-')
    options.channel = str2num([cell_list{1,i}(end-3) cell_list{1,i}(end-2)]);
  else
    options.channel = find(letter == cell_list{1,i}(end-1)); % convert letter to channel number
  end
  options.unit = str2num(cell_list{1,i}(end));
  
  %retrieve and analyze CLT
  cltinfo = findrecfromcell(cell_list{1,i}, '"CLT"');
  
%   % IF added 11 aug 2017 for getting rid of CLT elements in the list that were not sorted
%   if ~isempty(cltinfo) && ~any(cltinfo{end}.UnitsByElectrode{options.channel} == options.unit) % if unit is not present
%     cltinfo = []; % discard
%   end
  
  if ~isempty(cltinfo) && ~isempty(cltinfo{end})
    cltinfo = cltinfo{end};
   
    options.tag_masks = {'reference'};
    spikefilec = cltinfo.SpikeFile;
    paramsc    = cltinfo.exptparams;
    
    rrc = loadspikeraster(spikefilec, options);
    cltrates = paramsc.TrialObject.ReferenceHandle.ClickRate;
  else
    rrc=[];
  end
  
  
  
  if ~isempty(preidx) && length(actidx) == length(preidx)
    for j = 1 : length(actidx) % since num actives == passives we can do this
      options.tag_masks = {'SPECIAL-COLLAPSE-BOTH'};
      
      isol1 = cell_list{2,i}{preidx(j)}.SortInfo.Isolation;
      isol2 = cell_list{2,i}{actidx(j)}.SortInfo.Isolation;
      
      if nanmean(rrc(:)) > 2.5 % outliers! (VPr had some of them, I was filtering >10 before)
        
        rrc = [];
      end
      
      
      if (isol1 >= minisol || isol2 >= minisol) % && ~isempty(rrc) %11AUG17 removed this last condition
        spikefile1 = cell_list{2,i}{preidx(j)}.SpikeFile;
        spikefile2 = cell_list{2,i}{actidx(j)}.SpikeFile;
        params1    = cell_list{2,i}{preidx(j)}.exptparams;
        params2    = cell_list{2,i}{actidx(j)}.exptparams;
        
        r1 = loadspikeraster(spikefile1, options);
        r2 = loadspikeraster(spikefile2, options);
        
        %divide by delta t to have firing rate in spks/s
        r1 = (r1 .* psth_sr);
        r2 = (r2 .* psth_sr);
        if ~isempty(rrc), rrc=(rrc .* psth_sr); end
             
        
        % retrieve some parameters of current recording
        rrate    = params1.TrialObject.ReferenceHandle.ClickRate;
        trate    = params1.TrialObject.TargetHandle.ClickRate;
        clkonset = params1.TrialObject.ReferenceHandle.PreStimSilence + params1.TrialObject.ReferenceHandle.TorcDuration;
        clkoffset= respdur + clkonset;
        
        % determine if there was sound modulation & behavioral effect
        [h_temp h tmaxdiff mod passmod onsetidx] = reftar_stats_nonpar(r1,r2,params1, params2, psth_sr, alfa, 0);
                
        if h == 1 %&& mod == 1 % if modulated by sounds AND with significant behavioral effects %11AUG17 removed this condition
%           %% trim rasters to include only click part
%           r1 = r1((clkonset*psth_sr):(clkoffset*psth_sr), :, :);
%           r2 = r2((clkonset*psth_sr):(clkoffset*psth_sr), :, :);
          

          % save PSTHs
          rclk1(cnt,1:size(r1, 1)) = nanmean(squeeze(r1(:,:,1)), 2);
          tclk1(cnt,1:size(r1, 1)) = nanmean(squeeze(r1(:,:,2)), 2);
          rclk2(cnt,1:size(r2, 1)) = nanmean(squeeze(r2(:,:,1)), 2);
          tclk2(cnt,1:size(r2, 1)) = nanmean(squeeze(r2(:,:,2)), 2);
          
          if do_clt && ~isempty(rrc) 
            % find corresponding rates and psths in clt
            [~, rclosest] = min(abs(cltrates-rrate));
            [~, tclosest] = min(abs(cltrates-trate));
            
            rcltpsth = squeeze(rrc(:,:,rclosest)); %verificar que solo tiene 2 dim: tiempo x rep
            tcltpsth = squeeze(rrc(:,:,tclosest));
        
            
            % trim rasters
            cltonset= paramsc.TrialObject.ReferenceHandle.PreStimSilence;
            cltdur  = 1.2 + cltonset;
            
            %           rcltpsth = rcltpsth(cltonset*psth_sr:cltdur*psth_sr, :);
            %           tcltpsth = tcltpsth(cltonset*psth_sr:cltdur*psth_sr, :);
            
            rcltpsth = rcltpsth(1:cltdur*psth_sr, :);
            tcltpsth = tcltpsth(1:cltdur*psth_sr, :);
            
            % psths
            rcltpsth = nanmean(rcltpsth, 2); % verificar que sea de tamaño 1 x tiempo
            tcltpsth = nanmean(tcltpsth, 2);
            
            % save psths
            rclt(cnt, 1:length(rcltpsth)) = rcltpsth'; % verificar la forma de los psths
            tclt(cnt, 1:length(tcltpsth)) = tcltpsth';
            
            clt_cnt = clt_cnt + 1;
          end
          
          cellid_used{cnt} = cell_list{1,i}; % save cellID 
          parameters{cnt}  = params2; %save parameters
          %cellid_used{sel_cell_num+1} = cell_list{1,i}; % save cellID 
          sel_cell_num = sel_cell_num + 1; % count how many cells
        end
        
        cnt = cnt + 1;
      end
    end % for length(actidx)
    
  end % if preidx not empty and number of actives == passives
  
end % for ncells

% plot
t = [0, 1 : (respdur*psth_sr)]./ psth_sr; % had to add that +1 for dPEG, A1, dont understand why

% manually get rid of some outliers
% rclt([168, 171, 174], :) = NaN;
% tclt([168, 171, 174], :) = NaN;

% make matrices of a consistent size
rclk1 = rclk1(1:size(rclk1,1), 1:length(t));
tclk1 = tclk1(1:size(tclk1,1), 1:length(t));
rclk2 = rclk2(1:size(rclk2,1), 1:length(t));
tclk2 = tclk2(1:size(tclk2,1), 1:length(t));
rclt  = rclt (1:size(rclt, 1), 1:length(t));
tclt  = tclt (1:size(tclt, 1), 1:length(t));

% remove empty recordings
nonemptyrecs = find(~cellfun(@isempty, cellid_used));
rclk1 = rclk1(nonemptyrecs,:);
tclk1 = tclk1(nonemptyrecs,:);
rclk2 = rclk2(nonemptyrecs,:);
tclk2 = tclk2(nonemptyrecs,:);
rclt  = rclt(nonemptyrecs,:);
tclt  = tclt(nonemptyrecs,:);
cellid_used =cellid_used(nonemptyrecs);
parameters  = parameters(nonemptyrecs);

% % remove empty/NaN rows
% rclk1(~any(~isnan(rclk1),2),:)=[];
% tclk1(~any(~isnan(tclk1),2),:)=[];
% rclk2(~any(~isnan(rclk2),2),:)=[];
% tclk2(~any(~isnan(tclk2),2),:)=[];
% rclt (~any(~isnan(rclt),2),:) =[];
% tclt (~any(~isnan(tclt),2),:) =[];
% cellid_used = cellid_used(~cellfun('isempty', cellid_used));
% parameters = parameters(~cellfun('isempty', parameters));

% keyboard; %use only with A1 & dPEG data sets
% rclk1 = rclk1(:,1:109);
% rclk2 = rclk2(:,1:109);
% tclk1 = tclk1(:,1:109);
% tclk2 = tclk2(:,1:109);

mrclk1 = nanmean(rclk1);
mtclk1 = nanmean(tclk1);
mrclk2 = nanmean(rclk2);
mtclk2 = nanmean(tclk2);
mrclt  = nanmean(rclt);
mtclt  = nanmean(tclt);

[mrclk1e]= nanstd(rclk1) ./ sqrt(size(rclk1,1));
[mtclk1e]= nanstd(tclk1) ./ sqrt(size(tclk1,1));
[mrclk2e]= nanstd(rclk2) ./ sqrt(size(rclk2,1));
[mtclk2e]= nanstd(tclk2) ./ sqrt(size(tclk2,1));
[mrclte] = nanstd(rclt)  ./ sqrt(size(rclt,1));
[mtclte] = nanstd(tclt)  ./ sqrt(size(tclt,1));

% [a mrclk1e]= jackmeanerr(rclk1, 20);
% [a mtclk1e]= jackmeanerr(tclk1, 20);
% [a mrclk2e]= jackmeanerr(rclk2, 20);
% [a mtclk2e]= jackmeanerr(tclk2, 20);
% [a  mrclte] = jackmeanerr(rclt , 20);
% [a  mtclte] = jackmeanerr(tclt , 20);

% statistics

nbins = length(t);
nreps = size(rclk1, 1);

pclt = zeros(nbins, 1);
pclk1= zeros(nbins, 1);
pclk2= zeros(nbins, 1);

hclt = zeros(nbins, 1);
hclk1= zeros(nbins, 1);
hclk2= zeros(nbins, 1);

for i= 1: (3.2*psth_sr) % (3.2*psth_sr) respdur en vez de 3.2?
  %[hclt(i), pclt(i)] = ttest(rclt(:,i),tclt(:,i), alfa);
  [hclk1(i), pclk1(i)]= ttest(rclk1(:,i),tclk1(:,i), alfa);
  [hclk2(i), pclk2(i)]= ttest(rclk2(:,i),tclk2(:,i), alfa);
end

if do_clt && ~isempty(rclt) && ~isempty(tclt)
  for i=1: (cltdur*psth_sr)
    [hclt(i), pclt(i)] = ttest(rclt(:,i),tclt(:,i), alfa);
  end
    
  cltsig = find(pclt < alfa/nbins); % bonferroni correction
else
  cltsig = [];
end

clk1sig= find(pclk1 < alfa/nbins);
clk2sig= find(pclk2 < alfa/nbins);


figure;
set(gcf, 'Position', [190 323 1271 362]);
subplot(1,3,1);
if do_clt && ~isempty(rclt) && ~isempty(tclt)
  shadedErrorBar(t, mrclt, mrclte, 'b');
  hold on;
  shadedErrorBar(t, mtclt, mtclte, 'r');
  plot(t(cltsig), mtclt(cltsig), 'ro');
  hold off;
end

subplot(1,3,2);
shadedErrorBar(t, mrclk1, mrclk1e, 'b');
hold on;
shadedErrorBar(t, mtclk1, mtclk1e, 'r');
plot(t(clk1sig), mtclk1(clk1sig), 'ro');
hold off;

subplot(1,3,3);
shadedErrorBar(t, mrclk2, mrclk2e, 'b');
hold on;
shadedErrorBar(t, mtclk2, mtclk2e, 'r');
plot(t(clk2sig), mtclk2(clk2sig), 'ro');
hold off;

disp(['Number of cells: ', num2str(sel_cell_num), '/', num2str(cnt-1)]);

keyboard;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PFC/PFC_CLK.mat', ...
  'rclk1', 'rclk2', 'tclk1', 'tclk2', 'rclt', 'tclt','cellid_used', 't', 'parameters');

%% Organize and consolidate data from all areas into one file

clear all
load ('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/A1/A1_CLK.mat');

A1_cells_used = cellid_used;
A1_params     = parameters;
A1_Rpassive   = rclk1;
A1_Ractive    = rclk2;
A1_rCLT       = rclt;
A1_tCLT       = tclt;
A1_t          = t;
A1_Tpassive   = tclk1;
A1_Tactive    = tclk2;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/CLK_PSTHs.mat', 'A1*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/CLK_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/dPEG/dPEG_CLK.mat');
dPEG_cells_used = cellid_used;
dPEG_params     = parameters;
dPEG_Rpassive   = rclk1;
dPEG_Ractive    = rclk2;
dPEG_rCLT       = rclt;
dPEG_tCLT       = tclt;
dPEG_t          = t;
dPEG_Tpassive   = tclk1;
dPEG_Tactive    = tclk2;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/CLK_PSTHs.mat', 'A1*', 'dPEG*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/CLK_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/VPr/VPr_CLK.mat');
VPr_cells_used = cellid_used;
VPr_params     = parameters;
VPr_Rpassive   = rclk1;
VPr_Ractive    = rclk2;
VPr_rCLT       = rclt;
VPr_tCLT       = tclt;
VPr_t          = t;
VPr_Tpassive   = tclk1;
VPr_Tactive    = tclk2;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/CLK_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/CLK_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PFC/PFC_CLK.mat');
PFC_cells_used = cellid_used;
PFC_params     = parameters;
PFC_Rpassive   = rclk1;
PFC_Ractive    = rclk2;
PFC_rCLT       = rclt;
PFC_tCLT       = tclt;
PFC_t          = t;
PFC_Tpassive   = tclk1;
PFC_Tactive    = tclk2;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/CLK_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*', 'PFC*');
