% modified from SC_get_meanPSTH.m and ultimately from SC_meanPSTH_from_cell_list5a.m
% deleted all normalizations and plotting, just to be able to save the
% matrices with PSTHs and a list with cellids.
% It will enter debug mode once data has been collected.
% This version loads data from a naive animal (Tango), NO ACTIVES

% I'll try to make it work with tone and click data. But I'll leave the out
% of context clicks (CLT) for later...

% september 2017

doTORCpsth = 1;
doFTCpsth  = 1; clipFTC = 0; ftcmaxdur = 2;



accept_all_cells = 0;   %if 0 : check statistical significance of modulation from baseline and effects
%if 1 : consider all cells in cell_list

baphy_set_path;
dbopen;
dbstop if error;
%% load data -- output of ../database/getAreaRecordings
cd ('~/code/mycode/psth');
disp('loading cell list...');

area= 'VPr'; % A1, PEG, VPr or PFC
task= 'click'; % tone or click
normalize_sign = 0;
do_durfilt = 1; % bool: filter recordings by duration of target sound
% PFC: better use 0 (n=149) vs 1 (n=69)
% also svd_dPEG: 127 vs 77

if strcmp(task, 'click')
  doTORCpsth = 0;
  doFTCpsth  = 0;
end


switch area
  case 'VPr'
    if strcmp(task, 'tone')
      load ('Data/nPTD_VPr.mat');

    elseif strcmp(task, 'click')
      load ('Data/nCLK_VPr.mat');
    end
  case 'PFC'
    error('No naive data was collected in PFC');
  case 'PEG'
    if strcmp(task, 'tone')
      load ('Data/nPTD_PEG.mat');
    elseif strcmp(task, 'click')
      load ('Data/nCLK_PEG.mat');
    end
    
  case 'A1'
    if strcmp(task, 'tone')
      load ('Data/nPTD_A1.mat');
    elseif strcmp(task, 'click')
      load ('Data/nCLK_A1.mat');
    end
end


disp('Done Loading.');

% set duration of PSTHs and duration filter (for consistency)

durfilt    = 2;  %duration (s) of stimuli to use (do not consider silences here!)
presilence = 0.2;      % pre stimulus silence to include (s)
respdur    = 2.2;           % duration of response to consider (s) (from stim onset)

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

disp('Done retrieving cell isolation %');

%% find out how many passives and actives we have
passives=0;
actives=0;

for i = 1 : length(cell_list)
  passives = passives + sum(cell_list{3,i} == 'p');
  actives  = actives + sum(cell_list{3,i} == 'a');
end

%% add click discrminination direction info to the cell_list
if strcmp(task, 'click')
  row = cell(1,size(cell_list, 2));
  cell_list = [cell_list; row]; % create an empty row to store hl and lh discrim directions: u or d (going up or down, respectively)
  
  goingup = 0; % counters
  goingdown = 0;
  refrate =[]; tarrate = [];
  
  for i = 1 : length(cell_list)
    for j = 1 : length(cell_list{2,i})
      refrate = cell_list{2,i}{j}.exptparams.TrialObject.ReferenceHandle.ClickRate;
      tarrate = cell_list{2,i}{j}.exptparams.TrialObject.TargetHandle.ClickRate;
      if refrate < tarrate
        goingup = goingup + 1;
        cell_list{4,i} = [cell_list{4,i} 'u'];
      else
        goingdown = goingdown + 1;
        cell_list{4,i} = [cell_list{4,i} 'd'];
      end
    end
  end
  disp('Done determining click discrim directions.');
end


%% allocate PSTH matrices

cellid_used  = cell(1,passives);
parameters   = cell(1,passives);
passmods     = zeros(1,passives) .* NaN;
ftc_mods     = zeros(1,passives) .* NaN;
tor_mods     = zeros(1,passives) .* NaN;
sig_eff      = zeros(1,passives) .* NaN;

nbins = floor(options.rasterfs * (presilence+respdur));

RpassivePSTH = zeros(passives, nbins) .* NaN;
TpassivePSTH = zeros(passives, nbins) .* NaN;

if strcmp(task, 'click')
  RpassiveU = zeros(goingup, nbins) .* NaN;
  TpassiveU = zeros(goingup, nbins) .* NaN;
  RpassiveD = zeros(goingdown, nbins) .* NaN;
  TpassiveD = zeros(goingdown, nbins) .* NaN;
  cellid_usedU  = cell(1,goingup);
  parametersU   = cell(1,goingup);
  cellid_usedD  = cell(1,goingdown);
  parametersD   = cell(1,goingdown);  
  
end

if doTORCpsth
  TOR_PSTH     = zeros(ncell, nbins) .* NaN;
end

if doFTCpsth
  FTC_PSTH     = zeros(ncell, nbins) .* NaN;
end



%% retrieve and save PSTHs
passive_count = 0;
active_count  = 0;
tor_count     = 0;
ftc_count     = 0;

if strcmp(task, 'click')
  upcount = 0;
  downcount = 0;
end

%      %%%%%%%%%%%       MAIN LOOP - go through cell list  %%%%%%%%%%%%%%%      %
for i = 1 : length(cell_list)
  %cell_list{4,i}=[];
  
  % prepare some options
  if strcmp(cell_list{1,i}(end-1), '-')
    options.channel = str2num([cell_list{1,i}(end-3) cell_list{1,i}(end-2)]);
  else
    options.channel = find(letter == cell_list{1,i}(end-1)); % convert letter to channel number
  end
  options.unit = str2num(cell_list{1,i}(end));
  
  actidx = [];
  preidx = [];
  
  
  
  if length(cell_list{3,i}) > 1 % if more than one file was recorded
    
    preidx = find(cell_list{3,i} == 'p');
    

    for a = 1 : length(preidx)
   
      
     
      idx1       = preidx(a);
      spikefile1 = cell_list{2,i}{idx1}.SpikeFile;
      params1    = cell_list{2,i}{idx1}.exptparams;
      isol1      = cell_list{2,i}{idx1}.SortInfo.Isolation;
      if strcmp(task, 'click'), discdir = cell_list{4,i}(a); end
      
      if doFTCpsth
        taskfreq   = cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.Frequencies;
      end
        
      
                                                      
      if (isol1 >= minisolation) 
        
        if ~exist(spikefile1, 'file')
          keyboard;
          continue;
        end
        
        r1 = loadspikeraster(spikefile1, options);
        
        
        % go from spike count to firing rate
        r1 = r1 ./ (psth_sr/1000);
        
        
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
              r3 = loadspikeraster(TOR{1}.SpikeFile, TORoptions);
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
        
        if ~isempty(r1) 
          
          refpsth1 = r1(:,:,1);
          tarpsth1 = r1(:,:,2);
          

%           % test significance
%           [h_sig sig tmaxdiff mod passmod onsetidx] = reftar_stats_nonpar(r1, r2, params1, params2, (1000/psth_sr), alfa, 0);
          mod = modulated_from_baseline(r1, (1000/psth_sr), params1, alfa);

          if doFTCpsth && ~isempty(FTC)
            ftc_mod = modulated_from_baseline_singleref(sr4, FTCoptions.rasterfs, FTC{end}.exptparams, alfa); % check if response was modulated by sounds out of task context
          else
            ftc_mod = NaN;
          end
          
          if doTORCpsth && ~isempty(TOR)
            tor_mod = modulated_from_baseline_singleref(r3, TORoptions.rasterfs, TOR{1}.exptparams, alfa);
          else
            tor_mod = NaN;
          end
          % save results in cell_list
          %cell_list{5,i}.stats=struct('mod', mod, 'ftc_mod',ftc_mod, 'tor_mod', tor_mod);
          
    %      if strcmp(area, 'PFC')
            mod = 1; % override modulation testing 
%          end
          
          %% added to include sounds of just one duration
          if ~strcmp(params1.runclass, 'CLK')
            if do_durfilt && params1.TrialObject.TargetHandle.Duration == durfilt
              pass = 1;
            elseif ~do_durfilt
              pass = 1;
            else
              pass = 0;
            end
          
          else % if CLK
            
            if do_durfilt && params1.TrialObject.TargetHandle.TorcDuration + params1.TrialObject.TargetHandle.ClickDuration == durfilt
              pass = 1;
            elseif ~do_durfilt
              pass = 1;
            else
              pass = 0;
            end            
            
          end
          
          if pass ==1
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
            
            RpassivePSTH(passive_count + 1, 1:length(refpsth1)) = refpsth1;
            TpassivePSTH(passive_count + 1, 1:length(tarpsth1)) = tarpsth1;            
            
            if strcmp(task, 'click') && strcmp(discdir, 'u')
              upcount = upcount + 1;
              RpassiveU(upcount, 1:length(refpsth1)) = refpsth1;
              TpassiveU(upcount, 1:length(tarpsth1)) = tarpsth1;
              cellid_usedU{upcount} = cell_list{1,i};
              parametersU{upcount}  = params1;
              
            elseif strcmp(task, 'click') && strcmp(discdir, 'd')
              downcount = downcount + 1;  
              RpassiveD(downcount, 1:length(refpsth1)) = refpsth1;
              TpassiveD(downcount, 1:length(tarpsth1)) = tarpsth1;
              cellid_usedD{downcount} = cell_list{1,i};
              parametersD{downcount}  = params1;                          
            end
            
           
            passive_count = passive_count + 1;
            
            
            cellid_used{passive_count} = cell_list{1,i};
            parameters{passive_count}  = params1;
            ftc_mods(passive_count) = ftc_mod;
            tor_mods(passive_count) = tor_mod;

            
            if doTORCpsth && ~isempty(TOR) && ~isempty(r3)
              % TORC psths
              torpsth = nanmean(r3,2)';
              % trim psths to include only 0.2s silence and 1s sound response
              torsilence = TOR{1}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
              tordur = TOR{1}.exptparams.TrialObject.ReferenceHandle.Duration;
              
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
            
           
            
            cell_list{4,i} = [cell_list{4,i} 1]; % cell was selected
          else
            cell_list{4,i} = [cell_list{4,i} 0]; % cell was not selected
          end % if sig and mod == 1
        end % if r is empty
      end % if isolation
    end % for each active
    
  else % then there was only one file: must be active
    idx1 = 1;
    isol1 = cell_list{2,i}{1}.SortInfo.Isolation;
    if length(cell_list{3,i}) == 1 && strcmp(cell_list{3,i}, 'p') && isol1 >= minisolation %...
        % && strcmp(cell_list{2,i}{1}.IdentifierFull(1:3), 'was')
      spikefile1 = cell_list{2,i}{1}.SpikeFile;
      params1    = cell_list{2,i}{1}.exptparams;
      if strcmp(task, 'click'), discdir = cell_list{4,i}; end
      
      r1 = loadspikeraster(spikefile1, options);
      r1 = r1 ./ (psth_sr/1000);
      
      if doTORCpsth
        TORoptions = options;
        TORoptions.tag_masks= {'reference'};
        try
          TOR = findrecfromcell(cell_list{1,i}, '"TOR"');
          r3 = loadspikeraster(TOR{1}.SpikeFile, TORoptions);
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
        taskfreq   = cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.Frequencies;
        
        FTCoptions = options;
        FTCoptions.tag_masks= {'SPECIAL-TRIAL'};
        
        FTC = findrecfromcell(cell_list{1,i}, '"FTC"');
        if ~isempty(FTC)
          
          nftc = length(FTC);
          durftc = zeros(nftc,1);
          for nf = 1 : nftc
            durftc(nf) = FTC{nf}.exptparams.TrialObject.ReferenceHandle.Duration;
          end
          
          ftcidx = find(durftc > 0.1);
          
          try
            r4 = loadspikeraster(FTC{ftcidx(end)}.SpikeFile, FTCoptions);
            r4 = r4 ./ (psth_sr/1000);
          catch
            r4 = [];
            
          end
          
          
          % sort FTC raster by frequency
          nftctrials = size(r4,2);
          ftcfreqs = zeros(1,nftctrials); % this contains the frequency of each trial tone
          
          for f = 1: nftctrials
            [~, ~, Note] = evtimes(FTC{ftcidx(end)}.exptevents, 'Stim*', f);
            
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
            try
            r4 = sr4(:,closestf-2:closestf+2);
            r4freqs = sortedfreqs(closestf-2:closestf+2);
            catch
              r4 = sr4(:,1:5);
              r4freqs = sortedfreqs(1:5);
            end
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
        ftcsilence = FTC{ftcidx(end)}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
        ftcdur = FTC{ftcidx(end)}.exptparams.TrialObject.ReferenceHandle.Duration;
        
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
        FTC_PSTH(passive_count + 1, 1:length(ftcpsth)) = ftcpsth;
%         FTC_PSTH(ftc_count + 1, 1:length(ftcpsth)) = ftcpsth;
        ftc_count = ftc_count + 1;
      end
      
      refpsth1 = r1(:,:,1);
      tarpsth1 = r1(:,:,2);
      
      
      [mod, ~] = modulated_from_baseline(r1, (1000/psth_sr), params1, alfa);
      if doFTCpsth && ~isempty(FTC)
        ftc_mod = modulated_from_baseline_singleref(sr4, FTCoptions.rasterfs, FTC{ftcidx(end)}.exptparams, alfa); % check if response was modulated by sounds out of task context
      else
        ftc_mod = NaN;
      end
      
      if doTORCpsth && ~isempty(TOR)
        tor_mod = modulated_from_baseline_singleref(r3, TORoptions.rasterfs, TOR{1}.exptparams, alfa);
      else
        tor_mod = NaN;
      end
      
      if strcmp(params1.runclass, 'CLK')
        sounddur = params1.TrialObject.TargetHandle.TorcDuration + params1.TrialObject.TargetHandle.ClickDuration;
      else
        sounddur = params1.TrialObject.TargetHandle.Duration;
      end
      
      if do_durfilt && sounddur == durfilt
        pass = 1;
      elseif ~do_durfilt
        pass = 1;
      else
        pass = 0;
      end
      
      % override modulation testing:
      mod = 1;

      
      if (mod == 1 || accept_all_cells == 1) && pass == 1
        idx1 = 1;
        % passive
        refpsth1 = nanmean(refpsth1,2)';
        tarpsth1 = nanmean(tarpsth1,2)';
        % trim psths to include only 0.2s silence and 1s sound response
        refsilence = cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
        tarsilence = cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.PreStimSilence;
%         refdur = cell_list{2,i}{idx1}.exptparams.TrialObject.ReferenceHandle.Duration;
%         tardur = cell_list{2,i}{idx1}.exptparams.TrialObject.TargetHandle.Duration;
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
        
        rstartbin = refsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
        tstartbin = tarsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
        
        if rstartbin == 0, rstartbin =1; end
        if tstartbin == 0, tstartbin =1; end
        
        refpsth1 = refpsth1(rstartbin : (refsilence + respdur)*(1000/psth_sr));
        tarpsth1 = tarpsth1(tstartbin : (tarsilence + respdur)*(1000/psth_sr));
        
        if respdur_default > respdur % ugly fix for some shorter stimuli
          respdur = respdur_default;
        end
        
        RpassivePSTH(passive_count + 1, 1:length(refpsth1)) = refpsth1;
        TpassivePSTH(passive_count + 1, 1:length(tarpsth1)) = tarpsth1;
        cellid_used{passive_count + 1} = cell_list{1,i};
        parameters{passive_count+1}  = params1;
        ftc_mods(passive_count+1) = ftc_mod;
        tor_mods(passive_count+1) = tor_mod;
       
        if strcmp(task, 'click') && strcmp(discdir, 'u')
          upcount = upcount + 1;
          RpassiveU(upcount, 1:length(refpsth1)) = refpsth1;
          TpassiveU(upcount, 1:length(tarpsth1)) = tarpsth1;
          cellid_usedU{upcount} = cell_list{1,i};
          parametersU{upcount}  = params1;
          
        elseif strcmp(task, 'click') && strcmp(discdir, 'd')
          downcount = downcount + 1;
          RpassiveD(downcount, 1:length(refpsth1)) = refpsth1;
          TpassiveD(downcount, 1:length(tarpsth1)) = tarpsth1;
          cellid_usedD{downcount} = cell_list{1,i};
          parametersD{downcount}  = params1;
        end
                 
        
        if doTORCpsth && ~isempty(r3) && ~isempty(TOR) %  ~isempty(TOR{:})
          % TORC psths
          torpsth = nanmean(r3,2)';
          % trim psths to include only 0.2s silence and 1s sound response
          torsilence = TOR{1}.exptparams.TrialObject.ReferenceHandle.PreStimSilence;
          tordur = TOR{1}.exptparams.TrialObject.ReferenceHandle.Duration;
          
          torstartbin = torsilence*(1000/psth_sr) +1 - presilence*(1000/psth_sr);
          
          if torstartbin == 0, torstartbin =1; end
          
          torpsth = torpsth(torstartbin : (tordur)*(1000/psth_sr));
          TOR_PSTH(passive_count + 1, 1:length(torpsth)) = torpsth;
%           TOR_PSTH(tor_count + 1, 1:length(torpsth)) = torpsth;
          tor_count = tor_count + 1;
        end
        passive_count = passive_count + 1;
        
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
% DON'T remove all NaN rows, just trim it to used count rows
keyboard;

% trim matrices
RpassivePSTH = RpassivePSTH(1:passive_count, :);
TpassivePSTH = TpassivePSTH(1:passive_count, :);
if doFTCpsth, FTC_PSTH     = FTC_PSTH(1:passive_count, :); end
if doTORCpsth,TOR_PSTH     = TOR_PSTH(1:passive_count, 1:nbins); end
cellid_used  = cellid_used(1:passive_count);
parameters   = parameters(1:passive_count);
ftc_mods     = ftc_mods(1:passive_count);
tor_mods     = tor_mods(1:passive_count);

if strcmp(task, 'click')
  RpassiveU     = RpassiveU(1:upcount, :);
  TpassiveU     = TpassiveU(1:upcount, :);
  cellid_usedU  = cellid_usedU(1:upcount);
  parametersU   = parametersU(1:upcount);
  
  RpassiveD     = RpassiveD(1:downcount, :);
  TpassiveD     = TpassiveD(1:downcount, :);
  cellid_usedD  = cellid_usedD(1:downcount);
  parametersD   = parametersD(1:downcount); 
end

VPr_sr = options.rasterfs;

if strcmp(task, 'click')
  save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/A1/naive_VPr_CLK.mat', ...
    'VPr_sr','RpassivePSTH', 'TpassivePSTH', 'RpassiveU', 'TpassiveU', 'RpassiveD', 'TpassiveD',...
     'cellid_used','cellid_usedU', 'cellid_usedD', 't', 'parameters', 'parametersU', 'parametersD');  
else
  save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/A1/naive_A1_PTD.mat', ...
    'A1_sr','RpassivePSTH', 'TpassivePSTH', ...
    'FTC_PSTH', 'TOR_PSTH','cellid_used', 't', 'parameters', 'ftc_mods', 'tor_mods');
end


return;

%% make lists of click rates used

refrates = [];
tarrates = [];
tordur   = [];
clkdur   = [];

for i = 1 : length(parametersD)
  refrates = [refrates parametersD{i}.TrialObject.ReferenceHandle.ClickRate];
  tarrates = [tarrates parametersD{i}.TrialObject.TargetHandle.ClickRate];
  tordur   = [tordur parametersD{i}.TrialObject.ReferenceHandle.TorcDuration];
  clkdur   = [clkdur parametersD{i}.TrialObject.ReferenceHandle.ClickDuration];
end


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




figure;
subplot(3,3,[1 4]);
imagesc(FTC_PSTH, [0 90]);
%xlim([1 15]);
title('Tones alone');

subplot(3,3,7);
[TT TTe] = jackmeanerr(FTC_PSTH, 20);
shadedErrorBar(t, TT, TTe, 'r');
%xlim([-0.2 0.2]);
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

return;










%% ORGANIZING the individual areal PTD recordings

% 
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/A1/naive_A1_PTD.mat');
nA1_FTC = FTC_PSTH;
nA1_Rpassive = RpassivePSTH;
nA1_TOR = TOR_PSTH;
nA1_Tpassive = TpassivePSTH;
nA1_cells_used = cellid_used;
nA1_t = t;
nA1_sr = A1_sr;
nA1_params = parameters;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_TONE_PSTH.mat', 'nA1*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_TONE_PSTH.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/dPEG/naive_dPEG_PTD.mat');
ndPEG_FTC = FTC_PSTH;
ndPEG_Rpassive = RpassivePSTH;
ndPEG_TOR = TOR_PSTH;
ndPEG_Tpassive = TpassivePSTH;
ndPEG_cells_used = cellid_used;
ndPEG_t = t;
ndPEG_sr = dPEG_sr;
ndPEG_params = parameters;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_TONE_PSTH.mat', 'nA1*', 'ndPEG*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_TONE_PSTH.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/VPr/naive_VPr_PTD.mat');
nVPr_FTC = FTC_PSTH;
nVPr_Rpassive = RpassivePSTH;
nVPr_TOR = TOR_PSTH;
nVPr_Tpassive = TpassivePSTH;
nVPr_cells_used = cellid_used;
nVPr_t = t;
nVPr_sr = VPr_sr;
nVPr_params = parameters;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_TONE_PSTH.mat', 'nA1*', 'ndPEG*', 'nVPr*');

%% ORGANIZING the individual areal CLK recordings

% 
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/A1/naive_A1_CLK.mat');
nA1_Rpassive = RpassivePSTH;
nA1_Tpassive = TpassivePSTH;
nA1_cells_used = cellid_used;
nA1_RpassiveU = RpassiveU;
nA1_TpassiveU = TpassiveU;
nA1_cells_usedU = cellid_usedU;
nA1_RpassiveD = RpassiveD;
nA1_TpassiveD = TpassiveD;
nA1_cells_usedD = cellid_usedD;
nA1_t = t;
nA1_sr = A1_sr;
nA1_params = parameters;
nA1_paramsD = parametersD;
nA1_paramsU = parametersU;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_CLICK_PSTH.mat', 'nA1*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_CLICK_PSTH.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/dPEG/naive_dPEG_CLK.mat');
ndPEG_Rpassive = RpassivePSTH;
ndPEG_Tpassive = TpassivePSTH;
ndPEG_cells_used = cellid_used;
ndPEG_RpassiveU = RpassiveU;
ndPEG_TpassiveU = TpassiveU;
ndPEG_cells_usedU = cellid_usedU;
ndPEG_RpassiveD = RpassiveD;
ndPEG_TpassiveD = TpassiveD;
ndPEG_cells_usedD = cellid_usedD;
ndPEG_t = t;
ndPEG_sr = dPEG_sr;
ndPEG_params = parameters;
ndPEG_paramsD = parametersD;
ndPEG_paramsU = parametersU;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_CLICK_PSTH.mat', 'nA1*', 'ndPEG*');

clear all;
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_CLICK_PSTH.mat');
load('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/VPr/naive_VPr_CLK.mat');
nVPr_Rpassive = RpassivePSTH;
nVPr_Tpassive = TpassivePSTH;
nVPr_cells_used = cellid_used;
nVPr_RpassiveU = RpassiveU;
nVPr_TpassiveU = TpassiveU;
nVPr_cells_usedU = cellid_usedU;
nVPr_RpassiveD = RpassiveD;
nVPr_TpassiveD = TpassiveD;
nVPr_cells_usedD = cellid_usedD;
nVPr_t = t;
nVPr_sr = VPr_sr;
nVPr_params = parameters;
nVPr_paramsD = parametersD;
nVPr_paramsU = parametersU;
save('/home/delgueda/code/mycode/psth/Raw Data/30Hz-nomodfilt/naive_CLICK_PSTH.mat', 'nA1*', 'ndPEG*', 'nVPr*');

