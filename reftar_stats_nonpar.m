function [h_temp h tmaxdiff mod passmod onsetidx] = reftar_stats_nonpar(r1, r2, params1, params2, samplerate, alfa, debugplot)
% [h_temp h tmaxdiff mod onsetidx] = reftar_stats(r1, r2, params1, params2, samplerate, debugplot)
% will find significative differences between r1 and r2 (raster matrices,
% already smoothed and downsampled).
% h_temp: significance of changes in several time windows: onset (0-100ms),
% sustained(100-end), offset(end-end+0.1), task specific (end+.2 - end) and irrelevant sound (torc in CLK task).
% this function is still incomplete and works only in tone detection tasks

% necesito que retorne:
%    - retorne un ratio onset/sustained activity (medida de onset response, onsetidx)

%% definitions
sig_thresh = 2;

tonset = 0.1; % onset will by 0 to tonset (s). Sustained = tonset-end
toffset= 0.2; % offset resp = end-toffset

if ~exist('debugplot', 'var'), debugplot = 0; end

if strcmp(params1.runclass, 'PTD') || strcmp(params1.runclass, 'CCH')
  h_temp = [0,0,0,0; 0,0,0,0];  % 1st row:ref, 2nd row: tar
elseif strcmp(params1.runclass, 'CLK')
  h_temp = [0,0,0,0,0; 0,0,0,0,0];
end
    
h = 0;
mod = 0;
passmod = 0;
onsetidx = [];

%% check durations are the same
if ~strcmp(params1.runclass, 'CLK')
  if (params1.TrialObject.ReferenceHandle.PreStimSilence ~= params2.TrialObject.ReferenceHandle.PreStimSilence) ||...
      (params1.TrialObject.ReferenceHandle.Duration ~= params2.TrialObject.ReferenceHandle.Duration)
    warning('reftar_stats :  sound or silence durations are different between conditions!!');
    h_temp = []; h=0; tmaxdiff=0; mod=0; onsetidx=[];
    
    return;
  end
  
  
  if (params1.TrialObject.ReferenceHandle.PreStimSilence ~= params1.TrialObject.TargetHandle.PreStimSilence) ||...
      (params1.TrialObject.ReferenceHandle.Duration ~= params1.TrialObject.TargetHandle.Duration) || ...
      (params2.TrialObject.ReferenceHandle.PreStimSilence ~= params2.TrialObject.TargetHandle.PreStimSilence) ||...
      (params2.TrialObject.ReferenceHandle.Duration ~= params2.TrialObject.TargetHandle.Duration)
    warning('reftar_stats :  sound or silence durations are different between sound classes!!');
    h_temp = []; h=0; tmaxdiff=0; mod=0; onsetidx=[];
    
    return;
  end
  
  
  if (params1.TrialObject.ReferenceHandle.PostStimSilence ~= params2.TrialObject.ReferenceHandle.PostStimSilence) ||...
      (params1.TrialObject.TargetHandle.PostStimSilence ~= params2.TrialObject.TargetHandle.PostStimSilence)
    warning('reftar_stats :  poststim silence durations were different between conditions!!');
    
    minsize= min([size(r1,1), size(r2,1)]);
    r1 = r1(1:minsize, :, :);
    r2 = r2(1:minsize, :, :);
  end
  
else % if it is CLK
  
  if (params1.TrialObject.ReferenceHandle.PreStimSilence ~= params2.TrialObject.ReferenceHandle.PreStimSilence) ||...
      (params1.TrialObject.ReferenceHandle.TorcDuration ~= params2.TrialObject.ReferenceHandle.TorcDuration) || ...
      (params1.TrialObject.ReferenceHandle.ClickDuration ~= params2.TrialObject.ReferenceHandle.ClickDuration)
    warning('reftar_stats :  sound or silence durations are different between conditions!!');
    h_temp = []; h=0; tmaxdiff=0; mod=0; onsetidx=[];
    
    return;
  end
  
  
  if (params1.TrialObject.ReferenceHandle.PreStimSilence ~= params1.TrialObject.TargetHandle.PreStimSilence) ||...
      (params1.TrialObject.ReferenceHandle.TorcDuration ~= params1.TrialObject.TargetHandle.TorcDuration) || ...
      (params1.TrialObject.ReferenceHandle.ClickDuration ~= params1.TrialObject.TargetHandle.ClickDuration) || ...
      (params2.TrialObject.ReferenceHandle.PreStimSilence ~= params2.TrialObject.TargetHandle.PreStimSilence) ||...
      (params2.TrialObject.ReferenceHandle.TorcDuration ~= params2.TrialObject.TargetHandle.TorcDuration) || ...
      (params2.TrialObject.ReferenceHandle.ClickDuration ~= params2.TrialObject.TargetHandle.ClickDuration)
    warning('reftar_stats :  sound or silence durations are different between sound classes!!');
    h_temp = []; h=0; tmaxdiff=0; mod=0; onsetidx=[];
    
    return;
  end
  
  
  if (params1.TrialObject.ReferenceHandle.PostStimSilence ~= params2.TrialObject.ReferenceHandle.PostStimSilence) ||...
      (params1.TrialObject.TargetHandle.PostStimSilence ~= params2.TrialObject.TargetHandle.PostStimSilence)
    warning('reftar_stats :  poststim silence durations were different between conditions!!');
    
    minsize= min([size(r1,1), size(r2,1)]);
    r1 = r1(1:minsize, :, :);
    r2 = r2(1:minsize, :, :);
  end
  
end

%% check modulation from baseline
[h1 ~] = modulated_from_baseline(r1, samplerate, params1, alfa);
[h2 ~] = modulated_from_baseline(r2, samplerate, params2, alfa);

if h1==1 || h2==1
  mod = 1;
else
  mod = 0;
end

passmod = h1; % is there any modulation form baseline in the passive stimuli?

%% normalization (at least remove the baseline firing rate?)
r1 = remove_baseline_fr(r1, params1, samplerate);
r2 = remove_baseline_fr(r2, params2, samplerate);


%% assess significance across conditions (between r1 and r2)

% organize rasters
ref1 = r1(:,:,1)';
tar1 = r1(:,:,2)';

ref2 = r2(:,:,1)';
tar2 = r2(:,:,2)';



% remove NaN rows
ref1(~any(~isnan(ref1),2),:)=[];
tar1(~any(~isnan(tar1),2),:)=[];
ref2(~any(~isnan(ref2),2),:)=[];
tar2(~any(~isnan(tar2),2),:)=[];

% % remove NaN columns
% ref1 = ref1(:,all(~isnan(ref1)));
% ref2 = ref2(:,all(~isnan(ref2)));
% tar1 = tar1(:,all(~isnan(tar1)));
% tar2 = tar2(:,all(~isnan(tar2)));
% 
% if length(tar1) ~= length(tar2) || length(ref1) ~= length(ref2)
%   warning('reftar_stats :  durations between conditions are different, returning nothing');
%   h_temp = []; h=0; tmaxdiff=0; mod=0; onsetidx=[];
%   return;
% end

% significance between conditions
rnsmp = min([size(ref1,2) size(ref2,2)]);
minref = min([size(ref1,1) size(ref2, 1)]);

for ss = 1:rnsmp % don't try to do statistics in all nan rows (when different sound durations were used)
  if all(isnan(ref1(1:minref, ss)))
    rnsmp = ss-1;
    break;
  end
end


rp = zeros(1,rnsmp);
for rr = 1 : rnsmp;
  [rp(rr), ~]= signrank(ref1(1:minref ,rr), ref2(1:minref ,rr), 'alpha', alfa);
end

tnsmp = min([size(tar1,2) size(tar2,2)]);
mintar = min([size(tar1,1) size(tar2, 1)]);

for ss = 1:tnsmp % don't try to do statistics in all nan rows (when different sound durations were used)
  if all(isnan(tar1(1:mintar, ss)))
    tnsmp = ss-1;
    break;
  end
end


tp = zeros(1,tnsmp);
for tt = 1 : tnsmp;
  try
    [tp(tt), ~, t]= signrank(tar1(1:mintar ,tt), tar2(1:mintar ,tt), 'alpha', alfa);
  catch
    tp = ones(1,tnsmp); % if ranksum fails: discard data
  end
end

% bonferroni correction: fixed several errors 24 april 2017
% ref_sig = rp <= alfa/rnsmp; 
% tar_sig = tp <= alfa/tnsmp;
% 
% ref_sigidx = find(ref_sig);
% tar_sigidx = find(tar_sig);

% False Discovery Rate (Benjamini-Yekutieli)
rp_fdr = alfa * (1:rnsmp)/(rnsmp*sum(1./(1:rnsmp)));
k      = max(find(sort(rp) < rp_fdr));
if isempty(k)
  ref_sigidx = [];
  ref_sig = zeros(1,rnsmp);
else
  ref_sigidx = find(rp < rp_fdr(k));
  ref_sig    = rp < rp_fdr(k);
end

tp_fdr = alfa * (1:tnsmp)/(tnsmp*sum(1./(1:tnsmp)));
k      = max(find(sort(tp) < tp_fdr));
if isempty(k)
  tar_sigidx = [];
  tar_sig = zeros(1,tnsmp);
else
  tar_sigidx = find(tp < tp_fdr(k));
  tar_sig    = tp < tp_fdr(k);
end


% get mean (and sem) psths
[ref1 semref1] = jackmeanerr(ref1, 20);
[ref2 semref2] = jackmeanerr(ref2, 20);
[tar1 semtar1] = jackmeanerr(tar1, 20);
[tar2 semtar2] = jackmeanerr(tar2, 20);

midx = min([size(ref1, 2) size(ref2, 2)]);
dref = abs(ref2(1:midx) - ref1(1:midx));
midx = min([size(tar1, 2) size(tar2, 2)]);
dtar = abs(tar2(1:midx) - tar1(1:midx));

% % get mean (and sem) psths
% [ref1 semref1] = jackmeanerr(ref1, 20);
% [ref2 semref2] = jackmeanerr(ref2, 20);
% [tar1 semtar1] = jackmeanerr(tar1, 20);
% [tar2 semtar2] = jackmeanerr(tar2, 20);
% 
% % deltas (differences) will be active(2) - passive(1)
% dref = abs(ref2 - ref1);
% dtar = abs(tar2 - tar1);
% 
% % average SEM between conditions
% mse_ref = nanmean([semref1 ; semref2]);
% mse_tar = nanmean([semtar1 ; semtar2]);
% 
% ref_zscore = dref ./ mse_ref;
% tar_zscore = dtar ./ mse_tar;
% 
% ref_sig = ref_zscore >= sig_thresh;
% tar_sig = tar_zscore >= sig_thresh;
% 
% ref_sigidx = find(ref_zscore >= sig_thresh);
% tar_sigidx = find(tar_zscore >= sig_thresh);

t = ([1:length(ref2)] .* (1/samplerate)) - params1.TrialObject.ReferenceHandle.PreStimSilence; %assume all had same length



if max([max(abs(ref1)) max(abs(ref2)) max(abs(tar1)) max(abs(tar2))]) <= 10
  warning('reftar_stats : peak firing rate below 10 spikes/s, discarding!!');
  h_temp = []; h=0; tmaxdiff=0; mod=0; onsetidx=[];
  return;
  
end


%% time windows: 

%presilence and sustained should be different for tone of click tasks
%at the beginning check that silence and durations are the same for both
%sounds, then decide what to do

if strcmp(params1.runclass, 'CLK')
  rb0    = (params1.TrialObject.ReferenceHandle.PreStimSilence + ...
    params1.TrialObject.ReferenceHandle.TorcDuration) * samplerate; % bin where sound starts
  rb100  = rb0 + (0.1 * samplerate);      
  rboff  = rb0 + params1.TrialObject.ReferenceHandle.ClickDuration * samplerate;  % bin by sound offset
else
  rb0    = params1.TrialObject.ReferenceHandle.PreStimSilence * samplerate; % bin where sound starts
  rb100  = rb0 + (0.1 * samplerate);                                         % bin 100 ms after onset
  rboff  = rb0 + params1.TrialObject.ReferenceHandle.Duration * samplerate;  % bin by sound offset
end

rboff2 = rboff + (0.1 * samplerate);                                       % bin 100 ms after offset
%rb_end = rboff + params1.TrialObject.ReferenceHandle.PostStimSilence * samplerate;



if strcmp(params1.runclass, 'CLK')
  tb0    = (params1.TrialObject.TargetHandle.PreStimSilence+ params1.TrialObject.TargetHandle.TorcDuration) * samplerate; % bin where sound starts
  tb100  = tb0 + (0.1 * samplerate);                                         % bin 100 ms after onset  
  tboff  = tb0  + params1.TrialObject.TargetHandle.ClickDuration* samplerate;  % bin by sound offset
else
  tb0    = params1.TrialObject.TargetHandle.PreStimSilence * samplerate; % bin where sound starts
  tb100  = tb0 + (0.1 * samplerate);                                         % bin 100 ms after onset
  tboff  = tb0 + params1.TrialObject.TargetHandle.Duration * samplerate;  % bin by sound offset
end

tboff2 = tboff + (0.1 * samplerate);                                       % bin 100 ms after offset
%tb_end = tboff + params1.TrialObject.TargetHandle.PostStimSilence * samplerate;

% whole sound significance testing : more bins for alpha level and at least
% one series of 2 consecutive significant bins
if strcmp(params1.runclass, 'PTD') || strcmp(params1.runclass, 'CCH') || strcmp(params1.runclass, 'CLK')
  minsig = (rboff-rb0) * alfa;
  if sum(ref_sig(rb0:rboff)) > minsig || sum(tar_sig(tb0:tboff)) > minsig
    if ~isempty(findstr(ref_sig(rb0:rboff), [1,1])) || ~isempty(findstr(tar_sig(tb0:tboff), [1,1]))
      h = 1;
    else
      h = 0;
    end
  else
    h=0;
  end
  

  
end

% onset
minsig = 2
if sum(ref_sig(rb0:rb100)) >= minsig && ~isempty(findstr(ref_sig(rb0:rb100), [1,1]))
  h_temp(1,1) = 1;
else
  h_temp(1,1) = 0;
end

if sum(tar_sig(tb0:tb100)) >= minsig && ~isempty(findstr(tar_sig(rb0:rb100), [1,1]))
  h_temp(2,1) = 1;
else
  h_temp(2,1) = 0;
end  

% sustained

minsig = 2; % find at least 3 significant bins
if sum(ref_sig(rb100:rboff)) >= minsig && ~isempty(findstr(ref_sig(rb100:rboff), [1,1]))
    h_temp(1,2) = 1;
else
    h_temp(1,2) = 0;
end

if sum(tar_sig(tb100:tboff)) >= minsig && ~isempty(findstr(tar_sig(tb100:tboff), [1,1]))
    h_temp(2,2) = 1;
else
    h_temp(2,2) = 0;
end


% off response 
minsig = 2;
if sum(ref_sig(rboff:rboff2)) >= minsig && ~isempty(findstr(ref_sig(rboff:rboff2), [1,1]))
  h_temp(1,3) = 1;
else
  h_temp(1,3) = 0;
end

if sum(tar_sig(tboff:tboff2)) >= minsig && ~isempty(findstr(tar_sig(tboff:tboff2), [1,1]))
  h_temp(2,3) = 1;
else
  h_temp(2,3) = 0;
end


% after-off response (task specific?)
minsig = 2;
if sum(ref_sig(rboff2:end)) >= minsig && ~isempty(findstr(ref_sig(rboff2:end), [1,1]))
  h_temp(1,4) = 1;
else
  h_temp(1,4) = 0;
end

if sum(tar_sig(tboff2:end)) >= minsig && ~isempty(findstr(tar_sig(tboff2:end), [1,1]))
  h_temp(2,4) = 1;
else
  h_temp(2,4) = 0;
end


%% max change time bins [[add condition: only if significant!]]
if max(dtar) >= max(dref)
  tmaxdiff = t(dtar == max(dtar));
else 
  tmaxdiff = t(dref == max(dref));
end

%% onset / sustained ratio of effect
onsetidx = [];
% not sure how to do this yet...

tmaxdiff = tmaxdiff(end);


max([max(abs(ref1)) max(abs(ref2)) max(abs(tar1)) max(abs(tar2))])



%% debug use only : plot and wait for click
if debugplot == 1
  
  clf;
  refaxes = subplot(1,2,1);
  shadedErrorBar(t,ref1,semref1, '-b');
  hold on;
  shadedErrorBar(t, ref2, semref2, '-r');
  if ~isempty(ref_sigidx), plot(t(ref_sigidx), ref2(ref_sigidx), 'ob'); end
  hold off;
  
  subplot(1,2,2);
  shadedErrorBar(t, tar1, semtar1, '-b');
  hold on;
  shadedErrorBar(t, tar2, semtar2, '-r');
  if ~isempty(tar_sigidx), plot(t(tar_sigidx), tar2(tar_sigidx), 'ob'); end
  hold off;
  
  titlestr = sprintf('modulated: %i,   significant change: %i', mod, h);
  suptitle(titlestr);
  yy= ylim;
  set(refaxes, 'YLim', yy);
  h_temp
  tmaxdiff
  h
  mod
  k = waitforbuttonpress;
end
