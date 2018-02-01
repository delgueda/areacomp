% modified from SC_meanPSTH_from_cell_list5a.m
% deleted all normalizations and plotting, just to be able to save the
% matrices with PSTHs and a list with cellids.
% It will enter debug mode once data has been collected.


% April 2017

doTORCpsth = 1;
doFTCpsth  = 1; clipFTC = 1; ftcmaxdur = 0.1;

accept_all_cells = 0;   %if 0 : check statistical significance of modulation from baseline and effects
%if 1 : consider all cells in cell_list

baphy_set_path;
dbopen;
dbstop if error;
%% load data -- output of ../database/getAreaRecordings
cd ('~/code/mycode/psth');
disp('loading cell list...');

area= 'PEG';
task= 'tone'; % tone or click
normalize_sign = 0;
do_durfilt = 0; % bool: filter recordings by duration of target sound
% PFC: better use 0 (n=149) vs 1 (n=69)
% also svd_dPEG: 127 vs 77

switch area
  case 'VPr'
    if strcmp(task, 'tone')
      load ('Data/PTD-VPr.mat');
      cell_list1=cell_list;
      load ('Data/CCH-VPr.mat');
      cell_list=[cell_list1 cell_list];
    elseif strcmp(task, 'click')
      load ('Data/CLK-VPr.mat');
    end
  case 'PFC'
    if strcmp(task, 'tone')
      load ('Data/PTD-PFC.mat');
      cell_list1=cell_list;
      load ('Data/CCH-PFC.mat');
      cell_list=[cell_list1 cell_list];
    elseif strcmp(task, 'click')
      load ('Data/CLK-PFC.mat');
    end
  case 'PEG'
    if strcmp(task, 'tone')
      load ('Data/PTD-PEG.mat');
    elseif strcmp(task, 'click')
      load ('Data/CLK-PEG.mat');
    end
    
  case 'svdPEG'
    if strcmp(task, 'tone')
      load ('Data/dPEG_svd_cells.mat');
      accept_all_cells = 1;
    elseif strcmp(task, 'click')
      keyboard; % not ready
    end    
  case 'A1'
    if strcmp(task, 'tone')
      load ('Data/PTD-A1.mat');
    elseif strcmp(task, 'click')
      load ('Data/CLK-A1.mat');
    end
end
% load ('Data/svd_cell_list.mat');
% load ('Data/averagePSTH_VPr_cell_list.mat');

disp('Done.');

% set duration of PSTHs and duration filter (for consistency)
if strcmp(area, 'VPr')
  durfilt = 2;  %duration (s) of stimuli to use (do not consider silences here!)
  presilence = 0.2;      % pre stimulus silence to include (s)
  respdur = 3.2;           % duration of response to consider (s) (from stim onset)
else
  if strcmp(task, 'click')
    durfilt = 2;  %duration (s) of stimuli to use (do not consider silences here!)
    presilence = 0.2;      % pre stimulus silence to include (s)
    respdur = 3.2;           % duration of response to consider (s) (from stim onset)
  else
    durfilt = 1;  %duration (s) of stimuli to use (do not consider silences here!)
    presilence = 0.2;      % pre stimulus silence to include (s)
    respdur = 1.8;           % duration of response to consider (s) (from stim onset)
  end
end
respdur_default = respdur;




ncell = length(cell_list);

% smoothwin = 30; % sliding window (ms) era 30!
%psth_sr   = 50; % if 50 (ms) : downsample rasters to 20 Hz (PSTH resolution)


alfa = 0.05;

letter = ['abcdefghijklm'];
minisolation = 79.5;   % use only cells with isolation >= 79.5%


%% options for loadspikeraster
options.tag_masks = {'SPECIAL-COLLAPSE-BOTH'};
options.rasterfs = 30;  %20 Hz, 50 ms bins
options.includeprestim = 1;

psth_sr = 1000/ options.rasterfs; % bin length in ms

%% first: get isolation of cells

for i = 1 : length(cell_list)
  for j = 1 : length(cell_list{2, i})
    
    cell_list{2,i}{j}.SortInfo.Isolation = [];
    tmpquery = mysql(['SELECT isolation FROM gSingleRaw WHERE cellid="' cell_list{1,i} '" AND rawid="' num2str(cell_list{2,i}{j}.ID) '"']);
    if isempty(tmpquery) || isempty(tmpquery(1).isolation)
      tmpquery(1).isolation = 0; % assign zero isolation to those that have no data (discard)
    end
    
    cell_list{2,i}{j}.SortInfo.Isolation = tmpquery(1).isolation;
    
  end
end

%% find out how many passives and actives we have
passives=0;
actives=0;

for i = 1 : length(cell_list)
  passives = passives + sum(cell_list{3,i} == 'p');
  actives  = actives + sum(cell_list{3,i} == 'a');
end

%% allocate PSTH matrices

cellid_used  = cell(1,actives);
parameters   = cell(1,actives);
passmods     = zeros(1,actives) .* NaN;
ftc_mods     = zeros(1,actives) .* NaN;
tor_mods     = zeros(1,actives) .* NaN;
sig_eff      = zeros(1,actives) .* NaN;

nbins = options.rasterfs * (presilence+respdur);

RpassivePSTH = zeros(passives, nbins) .* NaN;
RactivePSTH  = zeros(actives, nbins) .* NaN;
TpassivePSTH = zeros(passives, nbins) .* NaN;
TactivePSTH  = zeros(actives, nbins) .* NaN;
if doTORCpsth
  TOR_PSTH     = zeros(ncell, nbins) .* NaN;
end

if doFTCpsth
  FTC_PSTH     = zeros(ncell, nbins) .* NaN;
end

%% allocate PSTH matrices - separated onset, sustained, offset and after-offset effects

RpassivePSTH1 = zeros(passives, nbins) .* NaN;
RactivePSTH1  = zeros(actives, nbins) .* NaN;
TpassivePSTH1 = zeros(passives, nbins) .* NaN;
TactivePSTH1  = zeros(actives, nbins) .* NaN;

RpassivePSTH2 = zeros(passives, nbins) .* NaN;
RactivePSTH2  = zeros(actives, nbins) .* NaN;
TpassivePSTH2 = zeros(passives, nbins) .* NaN;
TactivePSTH2  = zeros(actives, nbins) .* NaN;

RpassivePSTH3 = zeros(passives, nbins) .* NaN;
RactivePSTH3  = zeros(actives, nbins) .* NaN;
TpassivePSTH3 = zeros(passives, nbins) .* NaN;
TactivePSTH3  = zeros(actives, nbins) .* NaN;

RpassivePSTH4 = zeros(passives, nbins) .* NaN;
RactivePSTH4  = zeros(actives, nbins) .* NaN;
TpassivePSTH4 = zeros(passives, nbins) .* NaN;
TactivePSTH4  = zeros(actives, nbins) .* NaN;

%% retrieve and save PSTHs
passive_count = 0;
active_count  = 0;
tor_count     = 0;
ftc_count     = 0;

%      %%%%%%%%%%%       MAIN LOOP - go through cell list  %%%%%%%%%%%%%%%      %
for i = 1 : length(cell_list)
  cell_list{4,i}=[];
  
  % prepare some options
  if strcmp(cell_list{1,i}(end-1), '-')
    options.channel = str2num([cell_list{1,i}(end-3) cell_list{1,i}(end-2)]);
  else
    options.channel = find(letter == cell_list{1,i}(end-1)); % convert letter to channel number
  end
  options.unit = str2num(cell_list{1,i}(end));
  
  actidx = [];
  preidx = [];
  
  if length(cell_list{3,i}) > 1 % else : it must be a pfc active, don't test p-a behavior, just test modulation
    
    actidx = find(cell_list{3,i} == 'a');
    
    %for each active
    for a = 1: length(actidx)
      
      % find its passives ...
      if actidx(a) ~= 1 && cell_list{3,i}(actidx(a) - 1) == 'p'
        preidx = [preidx, actidx(a)-1];
      elseif actidx(a) ~= length(cell_list{3,i}) && cell_list{3,i}(actidx(a) + 1) == 'p'
        preidx = [preidx, actidx(a)+1];
      else
        preidx = find(cell_list{3,i} == 'p');
      end
      
    end % done selecting actives and passives done with this cell
    
    
    % if more than one file; for each active; [if isolation?]; if r not empty: significance; if significant; trim; save; count;
    % else (one file); check it's active & isolation; if r not empty; if modulation; trim, save, count
    
    % test significance here
    for a = 1 : length(actidx)
      
      idx2       = actidx(a);
      spikefile2 = cell_list{2,i}{idx2}.SpikeFile;
      params2    = cell_list{2,i}{idx2}.exptparams;
      isol2      = cell_list{2,i}{idx2}.SortInfo.Isolation;
      if doFTCpsth
        taskfreq   = cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.Frequencies;
      end
      
      if length(actidx) == length(preidx)
        idx1       = preidx(a);
        spikefile1 = cell_list{2,i}{idx1}.SpikeFile;
        params1    = cell_list{2,i}{idx1}.exptparams;
        isol1      = cell_list{2,i}{idx1}.SortInfo.Isolation;
        
      elseif length(preidx) == 1
        idx1       = preidx(1);
        spikefile1 = cell_list{2,i}{preidx(1)}.SpikeFile;
        params1    = cell_list{2,i}{preidx(1)}.exptparams;
        isol1      = cell_list{2,i}{preidx(1)}.SortInfo.Isolation;
        
      end
      
                                                          % quick way to
                                                          % discard data
                                                          % from lemon: && ~strcmp(cell_list{2,i}{1}.IdentifierFull(1:5), 'lemon')
      if (isol1 >= minisolation || isol2 >= minisolation) %&& strcmp(cell_list{2,i}{1}.IdentifierFull(1:3), 'was')
        
        if ~exist(spikefile1, 'file') || ~exist(spikefile2, 'file')
          keyboard;
          continue;
        end
        
        r1 = loadspikeraster(spikefile1, options);
        r2 = loadspikeraster(spikefile2, options);
        
        % go from spike count to firing rate
        r1 = r1 ./ (psth_sr/1000);
        r2 = r2 ./ (psth_sr/1000);
        
        if doTORCpsth
          TORoptions = options;
          TORoptions.tag_masks= {'reference'};
          
          TOR = findrecfromcell(cell_list{1,i}, '"TOR"');
          try
            if isempty(TOR{1})
              TOR=[];
            end
          end
          
          if ~isempty(TOR)
            try
              r3 = loadspikeraster(TOR{end}.SpikeFile, TORoptions);
              r3=reshape(r3, size(r3,1), (size(r3,2)*size(r3,3)));
              r3 = r3 ./ (psth_sr/1000);
            catch
              r3 = [];
            end
          else
            r3=[];
            
          end
        else
          r3=[];
        end
        
        if doFTCpsth
          FTCoptions = options;
          FTCoptions.tag_masks= {'SPECIAL-TRIAL'};
          
          FTC = findrecfromcell(cell_list{1,i}, '"FTC"');
          try
            if isempty(FTC{1})
              FTC=[];
            end
          end
          
          if ~isempty(FTC)
            try
              r4 = loadspikeraster(FTC{end}.SpikeFile, FTCoptions);
              r4 = r4 ./ (psth_sr/1000);
            catch
              r4 = [];
              
            end
            
            
            % sort FTC raster by frequency
            nftctrials = size(r4,2);
            ftcfreqs = zeros(1,nftctrials); % this contains the frequency of each trial tone
            
            for f = 1: nftctrials
              [~, ~, Note] = evtimes(FTC{end}.exptevents, 'Stim*', f);
              
              tmpnote = strsplit(Note{1}, ',');
              ftcfreqs(f) = str2num(tmpnote{2}); % extract the freq of each trial from exptevents
            end
            
            [sortedfreqs, sortidx] = sort(ftcfreqs);
            %srr4 = rr4(:, sortidx);
            sr4  = r4(:, sortidx); % sort raster by frequency presented
            
            % find FTC frequency closest to taskfreq
            % and restructure r4 with closest ftc freqs
            if length(taskfreq) == 1
              [~, closestf] = min(abs(sortedfreqs - taskfreq));
              if closestf <= 2
                r4 = sr4(:,1:5);
                r4freqs = sortedfreqs(1:5);
              else
                r4 = sr4(:,closestf-2:closestf+2);
                r4freqs = sortedfreqs(closestf-2:closestf+2);
              end
              % if there's an error: closestf might be by the upper limit
              % of sortedfreqs
            else
              closestf = [];
              r4 = []; r4freqs = [];
              for f = 1 : length(taskfreq)
                [~, tclosestf] = min(abs(sortedfreqs - taskfreq(f)));
                closestf = [closestf tclosestf];
                if tclosestf <= 2
                  r4 = [r4, sr4(:,1:5)];
                  r4freqs= [r4freqs, sortedfreqs(1:5)];
                else
                  r4 = [r4, sr4(:,tclosestf-2:tclosestf+2)];
                  r4freqs = [r4freqs, sortedfreqs(tclosestf-2:tclosestf+2)];
                end
                
              end
            end
            
            
          else % if no FTC was found
            r4=[];
            
          end
        end
        % r4 is the (already smoothed) raster (time x freqs) of the FTC freqencies closest to
        % task frequency (or frequencies)
        
        if ~isempty(r1) && ~isempty(r2)
          
          refpsth1 = r1(:,:,1);
          tarpsth1 = r1(:,:,2);
          
          refpsth2 = r2(:,:,1);
          tarpsth2 = r2(:,:,2);
          
          % test significance
          [h_sig sig tmaxdiff mod passmod onsetidx] = reftar_stats_nonpar(r1, r2, params1, params2, (1000/psth_sr), alfa, 0);
          
          if ~isempty(FTC)
            ftc_mod = modulated_from_baseline_singleref(sr4, FTCoptions.rasterfs, FTC{end}.exptparams, alfa); % check if response was modulated by sounds out of task context
          else
            ftc_mod = NaN;
          end
          
          if ~isempty(TOR)
            tor_mod = modulated_from_baseline_singleref(r3, TORoptions.rasterfs, TOR{end}.exptparams, alfa);
          else
            tor_mod = NaN;
          end
          % save results in cell_list
          cell_list{5,i}.stats=struct('mod', mod, 'signif', sig, 'ftc_mod',ftc_mod, 'tor_mod', tor_mod, 'passmod', passmod,'w_sig', h_sig, 'tmaxdiff', tmaxdiff);
          
    %      if strcmp(area, 'PFC')
            mod = 1; % override modulation testing for PFC!
%          end
          
          %% added to include sounds of just one duration
          if ~strcmp(params2.runclass, 'CLK')
            if do_durfilt && params2.TrialObject.TargetHandle.Duration == durfilt && ...
                params1.TrialObject.TargetHandle.Duration == durfilt
              pass = 1;
            elseif ~do_durfilt
              pass = 1;
            else
              pass = 0;
            end
          
          else % if CLK
            
            if do_durfilt && params2.TrialObject.TargetHandle.TorcDuration + params2.TrialObject.TargetHandle.ClickDuration == durfilt && ...
                params1.TrialObject.TargetHandle.TorcDuration + params1.TrialObject.TargetHandle.ClickDuration == durfilt
              pass = 1;
            elseif ~do_durfilt
              pass = 1;
            else
              pass = 0;
            end            
            
          end
          
          if ((sig == 1 && mod == 1) || accept_all_cells == 1) && pass ==1
            % passive
            refpsth1 = nanmean(refpsth1,2)';
            tarpsth1 = nanmean(tarpsth1,2)';
            % trim psths to include only 0.2s silence and 1s sound response
            refsilence = cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
            tarsilence = cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.PreStimSilence;
            if strcmp(params1.runclass, 'CLK')
              refdur = cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.TorcDuration + ...
                cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.ClickDuration + ...
                cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.PostStimSilence;
              tardur = cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.TorcDuration +...
                cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.ClickDuration + ...
                cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.PostStimSilence;              
            else
              refdur = cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.Duration + ...
                cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.PostStimSilence;
              tardur = cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.Duration + ...
                cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.PostStimSilence;
            end
            
            if respdur > refdur || respdur > tardur % ugly fix for some shorter stimuli
              respdur = min(refdur, tardur);
            end
            
            rstartbin = refsilence*(1000/psth_sr) + 1 - presilence*(1000/psth_sr);
            tstartbin = tarsilence*(1000/psth_sr) + 1 - presilence*(1000/psth_sr);
            
            %             if rstartbin == 0, rstartbin =1; end
            %             if tstartbin == 0, tstartbin =1; end
            
            refpsth1 = refpsth1(rstartbin : (refsilence + respdur)*(1000/psth_sr));
            tarpsth1 = tarpsth1(tstartbin : (tarsilence + respdur)*(1000/psth_sr));
            
            if respdur_default > respdur % ugly fix for some shorter stimuli
              respdur = respdur_default;
            end
            
            RpassivePSTH(active_count + 1, 1:length(refpsth1)) = refpsth1;
            TpassivePSTH(active_count + 1, 1:length(tarpsth1)) = tarpsth1;            
%             RpassivePSTH(passive_count + 1, 1:length(refpsth1)) = refpsth1;
%             TpassivePSTH(passive_count + 1, 1:length(tarpsth1)) = tarpsth1;
            
            if accept_all_cells == 0
              % save separate PSTHs by time window of effect
              if (h_sig(1,1) == 1 || h_sig(2,1) == 1) && (tmaxdiff >0 && tmaxdiff <= 0.1)
                RpassivePSTH1(passive_count + 1, 1:length(refpsth1)) = refpsth1;
                TpassivePSTH1(passive_count + 1, 1:length(tarpsth1)) = tarpsth1;
              end
            
            
              if (h_sig(1,2) == 1 || h_sig(2,2) == 1) && (tmaxdiff >0.1 && tmaxdiff <= 2)
                RpassivePSTH2(passive_count + 1, 1:length(refpsth1)) = refpsth1;
                TpassivePSTH2(passive_count + 1, 1:length(tarpsth1)) = tarpsth1;
              end
            
            
              if (h_sig(1,3) == 1 && h_sig(2,3) == 1) && (tmaxdiff >2 && tmaxdiff <= 2.1)
                RpassivePSTH3(passive_count + 1, 1:length(refpsth1)) = refpsth1;
                TpassivePSTH3(passive_count + 1, 1:length(tarpsth1)) = tarpsth1;
              end
            
            
              if (h_sig(1,4) == 1 && h_sig(2,4) == 1) && (tmaxdiff >2.1)
                RpassivePSTH4(passive_count + 1, 1:length(refpsth1)) = refpsth1;
                TpassivePSTH4(passive_count + 1, 1:length(tarpsth1)) = tarpsth1;
              end
            
    	    end
            passive_count = passive_count + 1;
            
            % active
            refpsth2 = nanmean(refpsth2,2)';
            tarpsth2 = nanmean(tarpsth2,2)';
            % trim psths to include only 0.2s silence and 1s sound response
            refsilence = cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
            tarsilence = cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.PreStimSilence;
            
            if strcmp(params2.runclass, 'CLK')
              refdur = cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.TorcDuration + ...
                cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.ClickDuration + ...
                cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.PostStimSilence;
              tardur = cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.TorcDuration + ...
                cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.ClickDuration + ...
                cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.PostStimSilence;
            else
              refdur = cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.Duration + ...
                cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.PostStimSilence;
              tardur = cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.Duration + ...
                cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.PostStimSilence;
            end
            
            if respdur > refdur || respdur > tardur % ugly fix for some shorter stimuli
              respdur = min(refdur, tardur);
            end
            
            rstartbin = refsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
            tstartbin = tarsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
            
            if rstartbin == 0, keyboard; end
            if tstartbin == 0, keyboard; end
            
            refpsth2 = refpsth2(rstartbin : (refsilence + respdur)*(1000/psth_sr));
            tarpsth2 = tarpsth2(tstartbin : (tarsilence + respdur)*(1000/psth_sr));
            
            if respdur_default > respdur % ugly fix for some shorter stimuli
              respdur = respdur_default;
            end
            
            RactivePSTH(active_count + 1, 1:length(refpsth2)) = refpsth2;
            TactivePSTH(active_count + 1, 1:length(tarpsth2)) = tarpsth2;
            
            if accept_all_cells == 0
            
              % save separate PSTHs by time window of effect
              if (h_sig(1,1) == 1 || h_sig(2,1) == 1) && (tmaxdiff >0 && tmaxdiff <= 0.1)
                RactivePSTH1(active_count + 1, 1:length(refpsth2)) = refpsth2;
                TactivePSTH1(active_count + 1, 1:length(tarpsth2)) = tarpsth2;
              end
            
            
              if (h_sig(1,2) == 1 || h_sig(2,2) == 1) && (tmaxdiff >0.1 && tmaxdiff <= 2)
                RactivePSTH2(active_count + 1, 1:length(refpsth2)) = refpsth2;
                TactivePSTH2(active_count + 1, 1:length(tarpsth2)) = tarpsth2;
              end
            
            
              if (h_sig(1,3) == 1 && h_sig(2,3) == 1) && (tmaxdiff >2 && tmaxdiff <= 2.1)
                RactivePSTH3(active_count + 1, 1:length(refpsth2)) = refpsth2;
                TactivePSTH3(active_count + 1, 1:length(tarpsth2)) = tarpsth2;
              end
            
            
              if (h_sig(1,4) == 1 && h_sig(2,4) == 1) && (tmaxdiff >2.1)
                RactivePSTH4(active_count + 1, 1:length(refpsth2)) = refpsth2;
                TactivePSTH4(active_count + 1, 1:length(tarpsth2)) = tarpsth2;
              end
            
            end
            cellid_used{active_count+1} = cell_list{1,i};
            parameters{active_count+1}  = params2;
            passmods(active_count+1) = passmod;
            ftc_mods(active_count+1) = ftc_mod;
            tor_mods(active_count+1) = tor_mod;
            sig_eff(active_count+1)  = sig;
            
%             active_count = active_count + 1; % moved to l. 588
            
            if doTORCpsth && ~isempty(TOR) && ~isempty(r3)
              % TORC psths
              torpsth = nanmean(r3,2)';
              % trim psths to include only 0.2s silence and 1s sound response
              torsilence = TOR{end}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
              tordur = TOR{end}.exptparams.TrialObject.ReferenceHandle.Duration;
              
              torstartbin = torsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
              
              if torstartbin == 0, torstartbin =1; end
              if tordur <= respdur
                dd = tordur;
              else
                dd = respdur;
              end
              torpsth = torpsth(torstartbin : (torsilence + dd)*(1000/psth_sr));
              TOR_PSTH(active_count + 1, 1:length(torpsth)) = torpsth;
%               TOR_PSTH(tor_count + 1, 1:length(torpsth)) = torpsth;
              tor_count = tor_count + 1;
            end
            
            
            
            if doFTCpsth && ~isempty(r4) % REPEAT THIS again after else (only 1 active)
              % FTC psths
              ftcpsth = nanmean(r4,2)';
              % trim psths to include only 0.2s silence and respdur sound response
              ftcsilence = FTC{end}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
              ftcdur = FTC{end}.exptparams.TrialObject.ReferenceHandle.Duration;
              
              ftcstartbin = ftcsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
              
              if ftcstartbin == 0, ftcstartbin =1; end
              
              if ftcstartbin < 0 % if silence recorded is less that what I want to save, pad with nans
                ftcpsth=[(zeros(1,abs(ftcstartbin)).* NaN), ftcpsth];
                ftcstartbin = 1;
              end
              
              if ftcdur <= respdur
                dd = ftcdur;
              else
                dd = respdur;
              end
              
              if clipFTC
                dd = ftcmaxdur;
              end
              
              ftcpsth = ftcpsth(ftcstartbin : (ftcsilence + dd)*(1000/psth_sr));
              FTC_PSTH(active_count + 1, 1:length(ftcpsth)) = ftcpsth;
%               FTC_PSTH(ftc_count + 1, 1:length(ftcpsth)) = ftcpsth;
              ftc_count = ftc_count + 1;
            end
            
            active_count = active_count + 1;
            
            cell_list{4,i} = [cell_list{4,i} 1]; % cell was selected
          else
            cell_list{4,i} = [cell_list{4,i} 0]; % cell was not selected
          end % if sig and mod == 1
        end % if r is empty
      end % if isolation
    end % for each active
    
  else % then there was only one file: must be active
    
    isol2 = cell_list{2,i}{1}.SortInfo.Isolation;
    if length(cell_list{3,i}) == 1 && strcmp(cell_list{3,i}, 'a') && isol2 >= minisolation %...
        % && strcmp(cell_list{2,i}{1}.IdentifierFull(1:3), 'was')
      spikefile2 = cell_list{2,i}{1}.SpikeFile;
      params2    = cell_list{2,i}{1}.exptparams;
      
      r2 = loadspikeraster(spikefile2, options);
      r2 = r2 ./ (psth_sr/1000);
      
      if doTORCpsth
        TORoptions = options;
        TORoptions.tag_masks= {'reference'};
        try
          TOR = findrecfromcell(cell_List{1,i}, '"TOR"');
          r3 = loadspikeraster(TOR{end}.SpikeFile, TORoptions);
          r3 = reshape(r3, size(r3,1), (size(r3,2)*size(r3,3)));
        catch
          r3=[];
        end
        if ~isempty(r3)
          r3 = r3 ./ (psth_sr/1000);
          
        end
        
      end
      
      % repeat the code for generating r4 here
      
      if doFTCpsth
        FTCoptions = options;
        FTCoptions.tag_masks= {'SPECIAL-TRIAL'};
        
        FTC = findrecfromcell(cell_list{1,i}, '"FTC"');
        if ~isempty(FTC)
          try
            r4 = loadspikeraster(FTC{end}.SpikeFile, FTCoptions);
            r4 = r4 ./ (psth_sr/1000);
          catch
            r4 = [];
            
          end
          
          
          % sort FTC raster by frequency
          nftctrials = size(r4,2);
          ftcfreqs = zeros(1,nftctrials); % this contains the frequency of each trial tone
          
          for f = 1: nftctrials
            [~, ~, Note] = evtimes(FTC{end}.exptevents, 'Stim*', f);
            
            tmpnote = strsplit(Note{1}, ',');
            ftcfreqs(f) = str2num(tmpnote{2}); % extract the freq of each trial from exptevents
          end
          
          [sortedfreqs, sortidx] = sort(ftcfreqs);
          %srr4 = rr4(:, sortidx);
          sr4  = r4(:, sortidx); % sort raster by frequency presented
          
          % find FTC frequency closest to taskfreq
          % and restructure r4 with closest ftc freqs
          if length(taskfreq) == 1
            [~, closestf] = min(abs(sortedfreqs - taskfreq));
            r4 = sr4(:,closestf-2:closestf+2);
            r4freqs = sortedfreqs(closestf-2:closestf+2);
            % IF ERROR: is taskfreq near the edges of sortedfreqs???
          else
            closestf = [];
            r4 = []; r4freqs = [];
            for f = 1 : length(taskfreq)
              [~, tclosestf] = min(abs(sortedfreqs - taskfreq(f)));
              closestf = [closestf tclosestf];
              r4 = [r4, sr4(:,tclosestf-2:tclosestf+2)];
              r4freqs = [r4freqs, sortedfreqs(tclosestf-2:tclosestf+2)];
              
            end
          end
          
          
        else % if no FTC was found
          r4=[];
          
        end
      end
      % r4 is the (already smoothed) raster (time x freqs) of the FTC freqencies closest to
      % task frequency (or frequencies)
      
      
      if doFTCpsth && ~isempty(r4) % REPEATING code before else
        % FTC psths
        ftcpsth = nanmean(r4,2)';
        % trim psths to include only 0.2s silence and respdur sound response
        ftcsilence = FTC{end}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
        ftcdur = FTC{end}.exptparams.TrialObject.ReferenceHandle.Duration;
        
        ftcstartbin = ftcsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
        
        if ftcstartbin == 0, ftcstartbin =1; end
        
        if ftcstartbin < 0 % if silence recorded is less that what I want to save, pad with nans
          ftcpsth=[(zeros(1,abs(ftcstartbin)).* NaN), ftcpsth];
          
        end
        
        if ftcdur <= respdur
          dd = ftcdur;
        else
          dd = respdur;
        end
        
        if clipFTC
          dd = ftcmaxdur;
        end
        
        ftcpsth = ftcpsth(ftcstartbin : (ftcsilence + dd)*(1000/psth_sr));
        FTC_PSTH(active_count + 1, 1:length(ftcpsth)) = ftcpsth;
%         FTC_PSTH(ftc_count + 1, 1:length(ftcpsth)) = ftcpsth;
        ftc_count = ftc_count + 1;
      end
      
      refpsth2 = r2(:,:,1);
      tarpsth2 = r2(:,:,2);
      
      
      [mod, ~] = modulated_from_baseline(r2, (1000/psth_sr), params2, alfa);
      if ~isempty(FTC)
        ftc_mod = modulated_from_baseline_singleref(sr4, FTCoptions.rasterfs, FTC{end}.exptparams, alfa); % check if response was modulated by sounds out of task context
      else
        ftc_mod = NaN;
      end
      
      if ~isempty(TOR)
        tor_mod = modulated_from_baseline_singleref(r3, TORoptions.rasterfs, TOR{end}.exptparams, alfa);
      else
        tor_mod = NaN;
      end
      
      if strcmp(params2.runclass, 'CLK')
        sounddur = params2.TrialObject.TargetHandle.TorcDuration + params2.TrialObject.TargetHandle.ClickDuration;
      else
        sounddur = params2.TrialObject.TargetHandle.Duration;
      end
      
      if do_durfilt && sounddur == durfilt
        pass = 1;
      elseif ~do_durfilt
        pass = 1;
      else
        pass = 0;
      end

      
      if (mod == 1 || accept_all_cells == 1) && pass == 1
        idx2 = 1;
        % active
        refpsth2 = nanmean(refpsth2,2)';
        tarpsth2 = nanmean(tarpsth2,2)';
        % trim psths to include only 0.2s silence and 1s sound response
        refsilence = cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
        tarsilence = cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.PreStimSilence;
        refdur = cell_list{2,i}{idx2}.exptparams.TrialObject.ReferenceHandle.Duration;
        tardur = cell_list{2,i}{idx2}.exptparams.TrialObject.TargetHandle.Duration;
        
        if respdur > refdur || respdur > tardur % ugly fix for some shorter stimuli
          respdur = min(refdur, tardur);
        end
        
        rstartbin = refsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
        tstartbin = tarsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
        
        if rstartbin == 0, rstartbin =1; end
        if tstartbin == 0, tstartbin =1; end
        
        refpsth2 = refpsth2(rstartbin : (refsilence + respdur)*(1000/psth_sr));
        tarpsth2 = tarpsth2(tstartbin : (tarsilence + respdur)*(1000/psth_sr));
        
        if respdur_default > respdur % ugly fix for some shorter stimuli
          respdur = respdur_default;
        end
        
        RactivePSTH(active_count + 1, 1:length(refpsth2)) = refpsth2;
        TactivePSTH(active_count + 1, 1:length(tarpsth2)) = tarpsth2;
        cellid_used{active_count + 1} = cell_list{1,i};
        parameters{active_count+1}  = params2;
        passmods(active_count+1) = passmod;
        ftc_mods(active_count+1) = ftc_mod;
        tor_mods(active_count+1) = tor_mod;
        sig_eff(active_count+1)  = sig;
%         active_count = active_count + 1; % moved to l. 792
        
        
        if doTORCpsth && ~isempty(r3) && ~isempty(TOR{:}) %  ~isempty(TOR{:})
          % TORC psths
          torpsth = nanmean(r3,2)';
          % trim psths to include only 0.2s silence and 1s sound response
          torsilence = TOR{end}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
          tordur = TOR{end}.exptparams.TrialObject.ReferenceHandle.Duration;
          
          torstartbin = torsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
          
          if torstartbin == 0, torstartbin =1; end
          
          torpsth = torpsth(torstartbin : (torsilence + respdur)*(1000/psth_sr));
          TOR_PSTH(active_count + 1, 1:length(torpsth)) = torpsth;
%           TOR_PSTH(tor_count + 1, 1:length(torpsth)) = torpsth;
          tor_count = tor_count + 1;
        end
        active_count = active_count + 1;
        
        cell_list{4,i} = [cell_list{4,i} 1];
      else
        cell_list{4,i} = [cell_list{4,i} 0];
      end
    end % if length(cell_list{3,i}) == 1 && strcmp(cell_list{3,i}, 'a') && isol2 >= minisolation
    
  end % if more than one file
  
  
end %for each cell : main loop
%% time vector
t0bin = 1 + (presilence*(1000/psth_sr)); % sound onset bin

if strcmp(area, 'VPr') || strcmp(area, 'PEG') || strcmp(area, 'PFC') || strcmp(area, 'A1')
  t = (([0:nbins-1] ./(1000/psth_sr)) - presilence) + (psth_sr/2)/1000;
  
else
  t = (([0:nbins] ./(1000/psth_sr)) - presilence) + (psth_sr/2)/1000;
end


%% stop here to save the PSTHs as mat files
% (PSTHs untouched, no normalization yet!)
% DON'T remove all NaN rows, just trim it to active_count rows
keyboard;

% trim matrices
RpassivePSTH = RpassivePSTH(1:active_count, :);
TpassivePSTH = TpassivePSTH(1:active_count, :);
RactivePSTH  = RactivePSTH (1:active_count, :);
TactivePSTH  = TactivePSTH (1:active_count, :);
FTC_PSTH     = FTC_PSTH(1:active_count, :);
TOR_PSTH     = TOR_PSTH(1:active_count, :);
cellid_used  = cellid_used(1:active_count);
parameters   = parameters(1:active_count);
passmods     = passmods(1:active_count);
ftc_mods     = ftc_mods(1:active_count);
tor_mods     = tor_mods(1:active_count);
sig_eff      = sig_eff(1:active_count);


dPEG_sr = options.rasterfs;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/dPEG/dPEG_PTD.mat', ...
  'dPEG_sr','RpassivePSTH', 'RactivePSTH', 'TpassivePSTH', 'TactivePSTH',...
  'FTC_PSTH', 'TOR_PSTH','cellid_used', 't', 'parameters', 'passmods', 'sig_eff', 'ftc_mods', 'tor_mods');

pcogcells = (sum(passmods == 0) / active_count) * 100;
disp([num2str(pcogcells) '% of neurons with effects did not respond to passive task sounds']);
return;
 % A1 9.4595   PEG 5.6701    svdPEG 13.3858     VPr 12.7753    PFC 76.5101
 
%% plot raw PSTHs and the mean PSTH

figure; % 
subplot(3,3,[1 4]);
imagesc(TOR_PSTH, [0 90]);
title('TORC alone');

subplot(3,3,7);
[RR RRe] = jackmeanerr(TOR_PSTH, 20);
shadedErrorBar(t, RR, RRe, 'b');
xlim([-0.2 3]);
ylim([5 20]);

subplot(3, 3, [2 5]);
imagesc(RpassivePSTH, [0 90]);
title('passive TORC in tone detect task');
subplot(3,3,8);
[RP RPe] = jackmeanerr(RpassivePSTH, 20);
shadedErrorBar(t, RP, RPe, 'b');
xlim([-0.2 3]);
ylim([5 20]);

subplot(3,3, [3 6]);
imagesc(RactivePSTH, [0 90]);
title('active TORC in tone detect task');
subplot(3,3, 9);
[RA RAe] = jackmeanerr(RactivePSTH, 20);
shadedErrorBar(t, RA, RAe, 'b');
xlim([-0.2 3]);
ylim([5 20]);




figure;
subplot(3,3,[1 4]);
imagesc(FTC_PSTH, [0 90]);
xlim([1 15]);
title('Tones alone');

subplot(3,3,7);
[TT TTe] = jackmeanerr(FTC_PSTH, 20);
shadedErrorBar(t, TT, TTe, 'r');
xlim([-0.2 0.2]);
ylim([5 20]);

subplot(3, 3, [2 5]);
imagesc(TpassivePSTH, [0 90]);
title('passive Tones in tone detect task');

subplot(3,3,8);
[TP TPe] = jackmeanerr(TpassivePSTH, 20);
shadedErrorBar(t, TP, TPe, 'r');
xlim([-0.2 3]);
ylim([5 20]);

subplot(3,3, [3 6]);
imagesc(TactivePSTH, [0 90]);
title('active Tones in tone detect task');
subplot(3,3, 9);
[TA TAe] = jackmeanerr(TactivePSTH, 20);
shadedErrorBar(t, TA, TAe, 'r');
xlim([-0.2 3]);
ylim([5 20]);

return;














%% CLK part
% taken from plot_click_contextcomp.m

cd ~/code/mycode/psth
dbopen;
baphy_set_path;

disp('loading cell_list... this might take several minutes');

load ('Data/CLK-PEG.mat');
disp('Done!');

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
    if isempty(tmpquery(1).isolation)
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
      
      
      if (isol1 >= minisol || isol2 >= minisol) && ~isempty(rrc)
        spikefile1 = cell_list{2,i}{preidx(j)}.SpikeFile;
        spikefile2 = cell_list{2,i}{actidx(j)}.SpikeFile;
        params1    = cell_list{2,i}{preidx(j)}.exptparams;
        params2    = cell_list{2,i}{actidx(j)}.exptparams;
        
        r1 = loadspikeraster(spikefile1, options);
        r2 = loadspikeraster(spikefile2, options);
        
        %divide by delta t to have firing rate in spks/s
        r1 = (r1 .* psth_sr);
        r2 = (r2 .* psth_sr);
        rrc=(rrc .* psth_sr);
        
        
        
        
        % retrieve some parameters of current recording
        rrate    = params1.TrialObject.ReferenceHandle.ClickRate;
        trate    = params1.TrialObject.TargetHandle.ClickRate;
        clkonset = params1.TrialObject.ReferenceHandle.PreStimSilence + params1.TrialObject.ReferenceHandle.TorcDuration;
        clkoffset= respdur + clkonset;
        
        % determine if there was sound modulation & behavioral effect
        [h_temp h tmaxdiff mod passmod onsetidx] = reftar_stats_nonpar(r1,r2,params1, params2, psth_sr, alfa, 0);
                
        if mod == 1 && h == 1 % if modulated by sounds AND with significant behavioral effects
%           %% trim rasters to include only click part
%           r1 = r1((clkonset*psth_sr):(clkoffset*psth_sr), :, :);
%           r2 = r2((clkonset*psth_sr):(clkoffset*psth_sr), :, :);
%           

          % save PSTHs
          rclk1(cnt,1:size(r1, 1)) = nanmean(squeeze(r1(:,:,1)), 2);
          tclk1(cnt,1:size(r1, 1)) = nanmean(squeeze(r1(:,:,2)), 2);
          rclk2(cnt,1:size(r2, 1)) = nanmean(squeeze(r2(:,:,1)), 2);
          tclk2(cnt,1:size(r2, 1)) = nanmean(squeeze(r2(:,:,2)), 2);
          
          % find corresponding rates and psths in clt
          [~, rclosest] = min(abs(cltrates-rrate));
          [~, tclosest] = min(abs(cltrates-trate));
          
          % smooth the rasters
          %       rcltpsth = conv2(ones(smoothwin,1)/smoothwin, 1, rrc(:,:,rclosest), 'same');
          %       tcltpsth = conv2(ones(smoothwin,1)/smoothwin, 1, rrc(:,:,tclosest), 'same');
          %       rcltpsth = rcltpsth(1:psth_sr:end,:);
          %       tcltpsth = tcltpsth(1:psth_sr:end,:);
          
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


% remove empty/NaN rows
rclk1(~any(~isnan(rclk1),2),:)=[];
tclk1(~any(~isnan(tclk1),2),:)=[];
rclk2(~any(~isnan(rclk2),2),:)=[];
tclk2(~any(~isnan(tclk2),2),:)=[];
rclt (~any(~isnan(rclt),2),:) =[];
tclt (~any(~isnan(tclt),2),:) =[];
cellid_used = cellid_used(~cellfun('isempty', cellid_used));
parameters = parameters(~cellfun('isempty', parameters));


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

for i= 1: (3.2*psth_sr)
  %[hclt(i), pclt(i)] = ttest(rclt(:,i),tclt(:,i), alfa);
  [hclk1(i), pclk1(i)]= ttest(rclk1(:,i),tclk1(:,i), alfa);
  [hclk2(i), pclk2(i)]= ttest(rclk2(:,i),tclk2(:,i), alfa);
end

for i=1: (cltdur*psth_sr)
  [hclt(i), pclt(i)] = ttest(rclt(:,i),tclt(:,i), alfa);
end


cltsig = find(pclt < alfa/nbins); % bonferroni correction
clk1sig= find(pclk1 < alfa/nbins);
clk2sig= find(pclk2 < alfa/nbins);


figure;
set(gcf, 'Position', [190 323 1271 362]);
subplot(1,3,1);
shadedErrorBar(t, mrclt, mrclte, 'b');
hold on;
shadedErrorBar(t, mtclt, mtclte, 'r');
plot(t(cltsig), mtclt(cltsig), 'ro');
hold off;

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
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz/dPEG/dPEG_CLK.mat', ...
  'rclk1', 'rclk2', 'tclk1', 'tclk2', 'rclt', 'tclt','cellid_used', 't', 'parameters');


%% ORGANIZING the individual areal PTD recordings
modcounts = 1;
% 
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/A1/A1_PTD.mat');
A1_FTC = FTC_PSTH;
A1_Ractive = RactivePSTH;
A1_Rpassive = RpassivePSTH;
A1_TOR = TOR_PSTH;
A1_Tactive = TactivePSTH;
A1_Tpassive = TpassivePSTH;
A1_cells_used = cellid_used;
A1_t = t;
A1_params = parameters;
if modcounts
  A1_pass_mods = passmods;
  A1_ftc_mods = ftc_mods;
  A1_tor_mods = tor_mods;
  A1_sig_eff  = sig_eff;
end
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat', 'A1*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/dPEG/dPEG_PTD.mat');
dPEG_FTC = FTC_PSTH;
dPEG_Ractive = RactivePSTH;
dPEG_Rpassive = RpassivePSTH;
dPEG_TOR = TOR_PSTH;
dPEG_Tactive = TactivePSTH;
dPEG_Tpassive = TpassivePSTH;
dPEG_cells_used = cellid_used;
dPEG_t = t;
dPEG_params = parameters;
if modcounts
  dPEG_pass_mods = passmods;
  dPEG_ftc_mods = ftc_mods;
  dPEG_tor_mods = tor_mods;
  dPEG_sig_eff  = sig_eff;
end
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat', 'A1*', 'dPEG*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/VPr/VPr_PTD.mat');
VPr_FTC = FTC_PSTH;
VPr_Ractive = RactivePSTH;
VPr_Rpassive = RpassivePSTH;
VPr_TOR = TOR_PSTH;
VPr_Tactive = TactivePSTH;
VPr_Tpassive = TpassivePSTH;
VPr_cells_used = cellid_used;
VPr_t = t;
VPr_params = parameters;
if modcounts
  VPr_pass_mods = passmods;
  VPr_ftc_mods = ftc_mods;
  VPr_tor_mods = tor_mods;
  VPr_sig_eff  = sig_eff;
end
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PFC/PFC_PTD.mat');
PFC_FTC = FTC_PSTH;
PFC_Ractive = RactivePSTH;
PFC_Rpassive = RpassivePSTH;
PFC_TOR = TOR_PSTH;
PFC_Tactive = TactivePSTH;
PFC_Tpassive = TpassivePSTH;
PFC_cells_used = cellid_used;
PFC_t = t;
PFC_params = parameters;
if modcounts
  PFC_pass_mods = passmods;
  PFC_ftc_mods = ftc_mods;
  PFC_tor_mods = tor_mods;
  PFC_sig_eff  = sig_eff;
end
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*', 'PFC*');


% add svd dPEG cells
clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/svd_dPEG/svd_dPEG_PTD.mat');
svd_dPEG_FTC = FTC_PSTH;
svd_dPEG_Ractive = RactivePSTH;
svd_dPEG_Rpassive = RpassivePSTH;
svd_dPEG_TOR = TOR_PSTH;
svd_dPEG_Tactive = TactivePSTH;
svd_dPEG_Tpassive = TpassivePSTH;
svd_dPEG_cells_used = cellid_used;
svd_dPEG_t = t;
svd_dPEG_params = parameters;
if modcounts
  svd_dPEG_pass_mods = passmods;
  svd_dPEG_ftc_mods = ftc_mods;
  svd_dPEG_tor_mods = tor_mods;
  svd_dPEG_sig_eff  = sig_eff;
end
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/PTD_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*', 'PFC*', 'svd_dPEG*');

%% add all dPEG (no filtering by statistics)

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/PTD_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/dPEG/alldPEG_PTD.mat');
alldPEG_FTC = FTC_PSTH;
alldPEG_Ractive = RactivePSTH;
alldPEG_Rpassive = RpassivePSTH;
alldPEG_TOR = TOR_PSTH;
alldPEG_Tactive = TactivePSTH;
alldPEG_Tpassive = TpassivePSTH;
alldPEG_cells_used = cellid_used;
alldPEG_t = t;
alldPEG_params = parameters;
save('/home/delgueda/code/mycode/psth/Raw Data/PTD_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*', 'PFC*', 'svd_dPEG*', 'alldPEG*');


%% compile 'all' data set (no behavioral effect)

clear all;
% load('/home/delgueda/code/mycode/psth/Raw Data/PTD_all_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz/A1/all_A1_PTD.mat');
A1_FTC = FTC_PSTH;
A1_Ractive = RactivePSTH;
A1_Rpassive = RpassivePSTH;
A1_TOR = TOR_PSTH;
A1_Tactive = TactivePSTH;
A1_Tpassive = TpassivePSTH;
A1_cells_used = cellid_used;
A1_t = t;
A1_params = parameters;
A1_sr = all_A1_sr;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PTD_all_PSTHs.mat', 'A1*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PTD_all_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz/dPEG/all_dPEG_PTD.mat');
dPEG_FTC = FTC_PSTH;
dPEG_Ractive = RactivePSTH;
dPEG_Rpassive = RpassivePSTH;
dPEG_TOR = TOR_PSTH;
dPEG_Tactive = TactivePSTH;
dPEG_Tpassive = TpassivePSTH;
dPEG_cells_used = cellid_used;
dPEG_t = t;
dPEG_params = parameters;
dPEG_sr = all_dPEG_sr;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PTD_all_PSTHs.mat', 'A1*', 'dPEG*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PTD_all_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz/VPr/all_VPr_PTD.mat');
VPr_FTC = FTC_PSTH;
VPr_Ractive = RactivePSTH;
VPr_Rpassive = RpassivePSTH;
VPr_TOR = TOR_PSTH;
VPr_Tactive = TactivePSTH;
VPr_Tpassive = TpassivePSTH;
VPr_cells_used = cellid_used;
VPr_t = t;
VPr_params = parameters;
VPr_sr = all_VPr_sr;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PTD_all_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*');


clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PTD_all_PSTHs.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PFC/all_PFC_PTD.mat');
PFC_FTC = FTC_PSTH;
PFC_Ractive = RactivePSTH;
PFC_Rpassive = RpassivePSTH;
PFC_TOR = TOR_PSTH;
PFC_Tactive = TactivePSTH;
PFC_Tpassive = TpassivePSTH;
PFC_cells_used = cellid_used;
PFC_t = t;
PFC_params = parameters;
PFC_sr = all_PFC_sr;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz/PTD_all_PSTHs.mat', 'A1*', 'dPEG*', 'VPr*', 'PFC*');
