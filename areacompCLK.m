% script to play with normalizations, statistics, and plotting of my data
% Click discrimination
% DEG 14-June-2017
% added naive september 2017

% l. 3519 : post-target vs behavior

% using data from SC_get_meanPSTH_CLK.m and
% naive_MU_context_analysis.m
% adapted from areacomp.m
runclass = 'CLK';
add_laminar = 1;
latstat = 0; % attentional effect latency test


A1p_bsl = [];
A1a_bsl = [];
dPEGp_bsl = [];
dPEGa_bsl = [];
VPrp_bsl = [];
VPra_bsl = [];
PFCp_bsl = [];
PFCa_bsl = [];

% load('./with_post_30Hz/CLK_PSTHs.mat');
load('./30Hz-nomodfilt/CLK_PSTHs.mat');

%% reassigning some dPEG ventral outliers to VPr
% these were recordings done by jonathan, which had features of VP in their
% tuning AND were recorded quite deep. Also, according to the measurements
% and photos, they seem to be the most ventral recordings assigned to dPEG,
% so given the depth and response features, they are more likely to be from
% VP than dPEG. (change done Dec 14, 2017)

% VPr_cells_used{284} = dPEG_cells_used{10};
% VPr_cells_used{285} = dPEG_cells_used{30};
% VPr_cells_used{286} = dPEG_cells_used{15};

dPEG_cells_used(10) = [];
dPEG_cells_used(30) = [];
dPEG_cells_used(15) = [];

% VPr_Ractive ([284,285,286],:) = dPEG_Ractive ([10, 30, 15],:);
% VPr_Rpassive([284,285,286],:) = dPEG_Rpassive([10, 30, 15],:);
% VPr_Tactive ([284,285,286],:) = dPEG_Tactive ([10, 30, 15],:);
% VPr_Tpassive([284,285,286],:) = dPEG_Tpassive([10, 30, 15],:);
dPEG_Ractive ([10, 30, 15],:) = [];
dPEG_Rpassive([10, 30, 15],:) = [];
dPEG_Tactive ([10, 30, 15],:) = [];
dPEG_Tpassive([10, 30, 15],:) = [];


% VPr_params{284} = dPEG_params{10};
% VPr_params{285} = dPEG_params{30};
% VPr_params{286} = dPEG_params{15};
dPEG_params(10) = [];
dPEG_params(30) = [];
dPEG_params(15) = [];

% VPr_rCLT([284,285,286],:) = dPEG_rCLT([10, 30, 15],:);
% VPr_tCLT([284,285,286],:) = dPEG_tCLT([10, 30, 15],:);
dPEG_rCLT([10, 30, 15],:) = [];
dPEG_tCLT([10, 30, 15],:) = [];


%%%%%% only 15
% % VPr_cells_used{284} = dPEG_cells_used{15};
% % VPr_Ractive (284,:) = dPEG_Ractive (15,:);
% % VPr_Rpassive(284,:) = dPEG_Rpassive(15,:);
% % VPr_Tactive (284,:) = dPEG_Tactive (15,:);
% % VPr_Tpassive(284,:) = dPEG_Tpassive(15,:);
% % VPr_params{284} = dPEG_params{15};
% % VPr_rCLT(284,:) = dPEG_rCLT(15,:);
% % VPr_tCLT(284,:) = dPEG_tCLT(15,:);
% % 
% % dPEG_cells_used(15) = [];
% % dPEG_Ractive (15,:) = [];
% % dPEG_Rpassive(15,:) = [];
% % dPEG_Tactive (15,:) = [];
% % dPEG_Tpassive(15,:) = [];
% % dPEG_params(15) = [];
% % dPEG_rCLT(15,:) = [];
% % dPEG_tCLT(15,:) = [];

%%


if add_laminar
  load('~/Data/laminar_sorted/A1klustacells.mat');
  A1_Rpassive  = [A1_Rpassive  ; rclk1];
  A1_Ractive   = [A1_Ractive   ; rclk2];
  A1_Tpassive  = [A1_Tpassive  ; tclk1];
  A1_Tactive   = [A1_Tactive   ; tclk2];
  A1_cells_used= [A1_cells_used; cellids_used];
  A1_params    = [A1_params    ; params];
end
  

add_naive = 1;




if add_naive
  load('./30Hz-nomodfilt/naive_CLICK_PSTH.mat');
  
  % get rid of 9 VPr cells with inconsistent torc/click durations
  nVPr_Rpassive   = nVPr_Rpassive(10:end,:);
  nVPr_Tpassive   = nVPr_Tpassive(10:end,:);
  nVPr_cells_used = nVPr_cells_used(10:end);
  nVPr_params     = nVPr_params(10:end);
  
  nVPr_RpassiveU   = nVPr_RpassiveU(4:end,:);
  nVPr_TpassiveU   = nVPr_TpassiveU(4:end,:);
  nVPr_cells_usedU = nVPr_cells_usedU(4:end);
  nVPr_paramsU     = nVPr_paramsU(4:end);
  
  nVPr_RpassiveD   = nVPr_RpassiveD(7:end,:);
  nVPr_TpassiveD   = nVPr_TpassiveD(7:end,:);
  nVPr_cells_usedD = nVPr_cells_usedD(7:end);
  nVPr_paramsD     = nVPr_paramsD(7:end);  
  
end

A1_sr  = 30;
dPEG_sr= 30;
VPr_sr = 30;
PFC_sr = 30;

% normalization settings
subtract_baseline = 1; % PAPER DEFAULT 1
normalize_each    = 0; % PAPER DEFAULT 0
indep_norm        = 0; % PAPER DEFAULT 0
popnorm           = 1; % PAPER DEFAULT 1

% selection rules
limit_single_tar  = 0; %1: only 1 tar, 0:any -1:exclude 1 tar
only_pos_tar      = 0; %1:only +, 0: any, -1:exclude + resp
pos_tar_act       = 1; %1:select positive targets resp in active, 0: passive target
svd_subset        = 0;

plot_ctx_figures  = 0;
scatwin = [60:80]; % bins to use for calculating the target-reference contrast in scatter and hist figures
nscatwin = [46:72];

%% make lists of parameters used and separate CLK responses based on disc direction

% A1
refrates = [];
tarrates = [];
A1presil = [];
A1tordur = [];
A1clkdur = [];

for i = 1 : length(A1_params)
  refrates = [refrates A1_params{i}.TrialObject.ReferenceHandle.ClickRate];
  tarrates = [tarrates A1_params{i}.TrialObject.TargetHandle.ClickRate];
  A1presil   = [A1presil A1_params{i}.TrialObject.ReferenceHandle.PreStimSilence];
  A1tordur   = [A1tordur A1_params{i}.TrialObject.ReferenceHandle.TorcDuration];
  A1clkdur   = [A1clkdur A1_params{i}.TrialObject.ReferenceHandle.ClickDuration];
end

A1_RpassiveU   = A1_Rpassive  (refrates - tarrates < 0,:);
A1_TpassiveU   = A1_Tpassive  (refrates - tarrates < 0,:);
A1_RactiveU    = A1_Ractive   (refrates - tarrates < 0,:);
A1_TactiveU    = A1_Tactive   (refrates - tarrates < 0,:);
A1_cells_usedU = A1_cells_used(refrates - tarrates < 0);

A1_RpassiveD   = A1_Rpassive  (refrates - tarrates > 0,:);
A1_TpassiveD   = A1_Tpassive  (refrates - tarrates > 0,:);
A1_RactiveD    = A1_Ractive   (refrates - tarrates > 0,:);
A1_TactiveD    = A1_Tactive   (refrates - tarrates > 0,:);
A1_cells_usedD = A1_cells_used(refrates - tarrates > 0);

% dPEG
refrates = [];
tarrates = [];
dPEGpresil = [];
dPEGtordur = [];
dPEGclkdur = [];

for i = 1 : length(dPEG_params)
  refrates = [refrates dPEG_params{i}.TrialObject.ReferenceHandle.ClickRate];
  tarrates = [tarrates dPEG_params{i}.TrialObject.TargetHandle.ClickRate];
  dPEGpresil   = [dPEGpresil dPEG_params{i}.TrialObject.ReferenceHandle.PreStimSilence];
  dPEGtordur   = [dPEGtordur dPEG_params{i}.TrialObject.ReferenceHandle.TorcDuration];
  dPEGclkdur   = [dPEGclkdur dPEG_params{i}.TrialObject.ReferenceHandle.ClickDuration];  
end

dPEG_RpassiveU   = dPEG_Rpassive  (refrates - tarrates < 0,:);
dPEG_TpassiveU   = dPEG_Tpassive  (refrates - tarrates < 0,:);
dPEG_RactiveU    = dPEG_Ractive   (refrates - tarrates < 0,:);
dPEG_TactiveU    = dPEG_Tactive   (refrates - tarrates < 0,:);
dPEG_cells_usedU = dPEG_cells_used(refrates - tarrates < 0);

dPEG_RpassiveD   = dPEG_Rpassive  (refrates - tarrates > 0,:);
dPEG_TpassiveD   = dPEG_Tpassive  (refrates - tarrates > 0,:);
dPEG_RactiveD    = dPEG_Ractive   (refrates - tarrates > 0,:);
dPEG_TactiveD    = dPEG_Tactive   (refrates - tarrates > 0,:);
dPEG_cells_usedD = dPEG_cells_used(refrates - tarrates > 0);

% VPr
refrates = [];
tarrates = [];
VPrpresil   = [];
VPrtordur   = [];
VPrclkdur   = [];

for i = 1 : length(VPr_params)
  refrates = [refrates VPr_params{i}.TrialObject.ReferenceHandle.ClickRate];
  tarrates = [tarrates VPr_params{i}.TrialObject.TargetHandle.ClickRate];
  VPrpresil   = [VPrpresil VPr_params{i}.TrialObject.ReferenceHandle.PreStimSilence];
  VPrtordur   = [VPrtordur VPr_params{i}.TrialObject.ReferenceHandle.TorcDuration];
  VPrclkdur   = [VPrclkdur VPr_params{i}.TrialObject.ReferenceHandle.ClickDuration];
end

VPr_RpassiveU   = VPr_Rpassive  (refrates - tarrates < 0,:);
VPr_TpassiveU   = VPr_Tpassive  (refrates - tarrates < 0,:);
VPr_RactiveU    = VPr_Ractive   (refrates - tarrates < 0,:);
VPr_TactiveU    = VPr_Tactive   (refrates - tarrates < 0,:);
VPr_cells_usedU = VPr_cells_used(refrates - tarrates < 0);

VPr_RpassiveD   = VPr_Rpassive  (refrates - tarrates > 0,:);
VPr_TpassiveD   = VPr_Tpassive  (refrates - tarrates > 0,:);
VPr_RactiveD    = VPr_Ractive   (refrates - tarrates > 0,:);
VPr_TactiveD    = VPr_Tactive   (refrates - tarrates > 0,:);
VPr_cells_usedD = VPr_cells_used(refrates - tarrates > 0);

% PFC
refrates = [];
tarrates = [];
PFCpresil= [];
PFCtordur= [];
PFCclkdur= [];

for i = 1 : length(PFC_params)
  refrates = [refrates PFC_params{i}.TrialObject.ReferenceHandle.ClickRate];
  tarrates = [tarrates PFC_params{i}.TrialObject.TargetHandle.ClickRate];
  PFCpresil   = [PFCpresil PFC_params{i}.TrialObject.ReferenceHandle.PreStimSilence];
  PFCtordur   = [PFCtordur PFC_params{i}.TrialObject.ReferenceHandle.TorcDuration];
  PFCclkdur   = [PFCclkdur PFC_params{i}.TrialObject.ReferenceHandle.ClickDuration];  
end

PFC_RpassiveU   = PFC_Rpassive  (refrates - tarrates < 0,:);
PFC_TpassiveU   = PFC_Tpassive  (refrates - tarrates < 0,:);
PFC_RactiveU    = PFC_Ractive   (refrates - tarrates < 0,:);
PFC_TactiveU    = PFC_Tactive   (refrates - tarrates < 0,:);
PFC_cells_usedU = PFC_cells_used(refrates - tarrates < 0);

PFC_RpassiveD   = PFC_Rpassive  (refrates - tarrates > 0,:);
PFC_TpassiveD   = PFC_Tpassive  (refrates - tarrates > 0,:);
PFC_RactiveD    = PFC_Ractive   (refrates - tarrates > 0,:);
PFC_TactiveD    = PFC_Tactive   (refrates - tarrates > 0,:);
PFC_cells_usedD = PFC_cells_used(refrates - tarrates > 0);

if add_naive
  nA1refrates = [];
  nA1tarrates = [];
  nA1presil = [];
  nA1tordur = [];
  nA1clkdur = [];
  
  for i = 1 : length(nA1_params)
    nA1refrates = [refrates nA1_params{i}.TrialObject.ReferenceHandle.ClickRate];
    nA1tarrates = [tarrates nA1_params{i}.TrialObject.TargetHandle.ClickRate];
    nA1presil   = [nA1presil nA1_params{i}.TrialObject.ReferenceHandle.PreStimSilence];
    nA1tordur   = [nA1tordur nA1_params{i}.TrialObject.ReferenceHandle.TorcDuration];
    nA1clkdur   = [nA1clkdur nA1_params{i}.TrialObject.ReferenceHandle.ClickDuration];
  end

  ndPEGrefrates = [];
  ndPEGtarrates = [];
  ndPEGpresil = [];
  ndPEGtordur = [];
  ndPEGclkdur = [];
  
  for i = 1 : length(ndPEG_params)
    ndPEGrefrates = [refrates ndPEG_params{i}.TrialObject.ReferenceHandle.ClickRate];
    ndPEGtarrates = [tarrates ndPEG_params{i}.TrialObject.TargetHandle.ClickRate];
    ndPEGpresil   = [ndPEGpresil ndPEG_params{i}.TrialObject.ReferenceHandle.PreStimSilence];
    ndPEGtordur   = [ndPEGtordur ndPEG_params{i}.TrialObject.ReferenceHandle.TorcDuration];
    ndPEGclkdur   = [ndPEGclkdur ndPEG_params{i}.TrialObject.ReferenceHandle.ClickDuration];
  end  
  
  nVPrrefrates = [];
  nVPrtarrates = [];
  nVPrpresil = [];
  nVPrtordur = [];
  nVPrclkdur = [];
  
  for i = 1 : length(nVPr_params)
    nVPrrefrates = [refrates nVPr_params{i}.TrialObject.ReferenceHandle.ClickRate];
    nVPrtarrates = [tarrates nVPr_params{i}.TrialObject.TargetHandle.ClickRate];
    nVPrpresil   = [nVPrpresil nVPr_params{i}.TrialObject.ReferenceHandle.PreStimSilence];
    nVPrtordur   = [nVPrtordur nVPr_params{i}.TrialObject.ReferenceHandle.TorcDuration];
    nVPrclkdur   = [nVPrclkdur nVPr_params{i}.TrialObject.ReferenceHandle.ClickDuration];
  end
  
  
end

% NORMALIZATION
%% A1

nbins = size(A1_Rpassive, 2);

for i=1:size(A1_Rpassive, 1)
  % baseline mean removal
  
  tonset = 0.2 * A1_sr; % presilence was fixed at 0.2s
  
  A1p_bsl = [A1p_bsl nanmean([A1_Rpassive(i,1:tonset) A1_Tpassive(i,1:tonset)])];
  A1a_bsl = [A1a_bsl nanmean([A1_Ractive(i,1:tonset)  A1_Tactive(i,1:tonset)])];  
  
  if subtract_baseline
    bsl = A1_Rpassive(i,1:tonset);
    A1_Rpassive(i,:) = A1_Rpassive(i,:) - nanmean(bsl);
    
    bsl = A1_Tpassive(i,1:tonset);
    A1_Tpassive(i,:) = A1_Tpassive(i,:) - nanmean(bsl);
    
    bsl = A1_Ractive(i,1:tonset);
    A1_Ractive(i,:) = A1_Ractive(i,:) - nanmean(bsl);
    
    bsl = A1_Tactive(i,1:tonset);
    A1_Tactive(i,:) = A1_Tactive(i,:) - nanmean(bsl);
    
%     bsl = A1_rCLT(i,1:tonset);
%     A1_rCLT(i,:) = A1_rCLT(i,:) - nanmean(bsl);
%     
%     bsl = A1_tCLT(i,1:tonset);
%     A1_tCLT(i,:) = A1_tCLT(i,:) - nanmean(bsl);
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([A1_Rpassive(i,:) A1_Ractive(i,:) A1_Tpassive(i,:) A1_Tactive(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([A1_Rpassive(i,:) A1_Tpassive(i,:)]);
      fmaxA = max([A1_Ractive(i,:) A1_Tactive(i,:)]);
    end
    
    A1_Rpassive(i,:) = A1_Rpassive(i,:) ./ fmaxP;
    A1_Tpassive(i,:) = A1_Tpassive(i,:) ./ fmaxP;
    A1_Ractive(i,:) = A1_Ractive(i,:)   ./ fmaxA;
    A1_Tactive(i,:) = A1_Tactive(i,:)   ./ fmaxA;
    
%     A1_rCLT(i,:)     = A1_rCLT(i,:)       ./ fmaxP;
%     A1_tCLT(i,:)     = A1_tCLT(i,:)       ./ fmaxP;


  end
    
end

% A1 UP
for i=1:size(A1_RpassiveU, 1)
  % baseline mean removal
  
  tonset = 0.2 * A1_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = A1_RpassiveU(i,1:tonset);
    A1_RpassiveU(i,:) = A1_RpassiveU(i,:) - nanmean(bsl);
    
    bsl = A1_TpassiveU(i,1:tonset);
    A1_TpassiveU(i,:) = A1_TpassiveU(i,:) - nanmean(bsl);
    
    bsl = A1_RactiveU(i,1:tonset);
    A1_RactiveU(i,:) = A1_RactiveU(i,:) - nanmean(bsl);
    
    bsl = A1_TactiveU(i,1:tonset);
    A1_TactiveU(i,:) = A1_TactiveU(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([A1_RpassiveU(i,:) A1_RactiveU(i,:) A1_TpassiveU(i,:) A1_TactiveU(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([A1_RpassiveU(i,:) A1_TpassiveU(i,:)]);
      fmaxA = max([A1_RactiveU(i,:) A1_TactiveU(i,:)]);
    end
    
    A1_RpassiveU(i,:) = A1_RpassiveU(i,:) ./ fmaxP;
    A1_TpassiveU(i,:) = A1_TpassiveU(i,:) ./ fmaxP;
    A1_RactiveU(i,:) = A1_RactiveU(i,:)   ./ fmaxA;
    A1_TactiveU(i,:) = A1_TactiveU(i,:)   ./ fmaxA;
    
  end
    
end

% A1 DOWN
for i=1:size(A1_RpassiveD, 1)
  % baseline mean removal
  
  tonset = 0.2 * A1_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = A1_RpassiveD(i,1:tonset);
    A1_RpassiveD(i,:) = A1_RpassiveD(i,:) - nanmean(bsl);
    
    bsl = A1_TpassiveD(i,1:tonset);
    A1_TpassiveD(i,:) = A1_TpassiveD(i,:) - nanmean(bsl);
    
    bsl = A1_RactiveD(i,1:tonset);
    A1_RactiveD(i,:) = A1_RactiveD(i,:) - nanmean(bsl);
    
    bsl = A1_TactiveD(i,1:tonset);
    A1_TactiveD(i,:) = A1_TactiveD(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([A1_RpassiveD(i,:) A1_RactiveD(i,:) A1_TpassiveD(i,:) A1_TactiveD(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([A1_RpassiveD(i,:) A1_TpassiveD(i,:)]);
      fmaxA = max([A1_RactiveD(i,:) A1_TactiveD(i,:)]);
    end
    
    A1_RpassiveD(i,:) = A1_RpassiveD(i,:) ./ fmaxP;
    A1_TpassiveD(i,:) = A1_TpassiveD(i,:) ./ fmaxP;
    A1_RactiveD(i,:) = A1_RactiveD(i,:)   ./ fmaxA;
    A1_TactiveD(i,:) = A1_TactiveD(i,:)   ./ fmaxA;
    
  end
    
end

if popnorm
  [A1_Rpassive, A1_Ractive, A1_Tpassive, A1_Tactive] = popmaxnorm(A1_Rpassive, A1_Ractive, A1_Tpassive, A1_Tactive);
  [A1_RpassiveU, A1_RactiveU, A1_TpassiveU, A1_TactiveU] = popmaxnorm(A1_RpassiveU, A1_RactiveU, A1_TpassiveU, A1_TactiveU);
  [A1_RpassiveD, A1_RactiveD, A1_TpassiveD, A1_TactiveD] = popmaxnorm(A1_RpassiveD, A1_RactiveD, A1_TpassiveD, A1_TactiveD);
end


if add_naive
  % naive (tango) data
  tonset = 0.2 * nA1_sr; % presilence was fixed at 0.2s
  
  for i=1:size(nA1_Rpassive, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = nA1_Rpassive(i,1:tonset);
      nA1_Rpassive(i,:) = nA1_Rpassive(i,:) - nanmean(bsl);
      
      bsl = nA1_Tpassive(i,1:tonset);
      nA1_Tpassive(i,:) = nA1_Tpassive(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([nA1_Rpassive(i,:) nA1_Tpassive(i,:)]);
      nA1_Rpassive(i,:) = nA1_Rpassive(i,:) ./ fmax;
      nA1_Tpassive(i,:) = nA1_Tpassive(i,:) ./ fmax;
    end
  end
  
  % click discrim going UP
  for i=1:size(nA1_RpassiveU, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = nA1_RpassiveU(i,1:tonset);
      nA1_RpassiveU(i,:) = nA1_RpassiveU(i,:) - nanmean(bsl);
      
      bsl = nA1_TpassiveU(i,1:tonset);
      nA1_TpassiveU(i,:) = nA1_TpassiveU(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([nA1_RpassiveU(i,:) nA1_TpassiveU(i,:)]);
      nA1_RpassiveU(i,:) = nA1_RpassiveU(i,:) ./ fmax;
      nA1_TpassiveU(i,:) = nA1_TpassiveU(i,:) ./ fmax;
    end
  end
  
  % click discrim going DOWN
  for i=1:size(nA1_RpassiveD, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = nA1_RpassiveD(i,1:tonset);
      nA1_RpassiveD(i,:) = nA1_RpassiveD(i,:) - nanmean(bsl);
      
      bsl = nA1_TpassiveD(i,1:tonset);
      nA1_TpassiveD(i,:) = nA1_TpassiveD(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([nA1_RpassiveD(i,:) nA1_TpassiveD(i,:)]);
      nA1_RpassiveD(i,:) = nA1_RpassiveD(i,:) ./ fmax;
      nA1_TpassiveD(i,:) = nA1_TpassiveD(i,:) ./ fmax;
    end
  end  
  
  if popnorm
    [nA1_Rpassive, nA1_Tpassive] = popmaxnorm(nA1_Rpassive, nA1_Tpassive);
    [nA1_RpassiveU, nA1_TpassiveU] = popmaxnorm(nA1_RpassiveU, nA1_TpassiveU);
    [nA1_RpassiveD, nA1_TpassiveD] = popmaxnorm(nA1_RpassiveD, nA1_TpassiveD);
  end
end

%% dPEG

nbins = size(dPEG_Rpassive, 2);

for i=1:size(dPEG_Rpassive, 1)
  % baseline mean removal
  
  tonset = 0.2 * dPEG_sr; % presilence was fixed at 0.2s
  %bsl = zeros(1,nbins).*NaN;
  
  dPEGp_bsl = [dPEGp_bsl nanmean([dPEG_Rpassive(i,1:tonset) dPEG_Tpassive(i,1:tonset)])];
  dPEGa_bsl = [dPEGa_bsl nanmean([dPEG_Ractive(i,1:tonset)  dPEG_Tactive(i,1:tonset)])];  
  
  if subtract_baseline
    bsl = dPEG_Rpassive(i,1:tonset);
    dPEG_Rpassive(i,:) = dPEG_Rpassive(i,:) - nanmean(bsl);
    
    bsl = dPEG_Tpassive(i,1:tonset);
    dPEG_Tpassive(i,:) = dPEG_Tpassive(i,:) - nanmean(bsl);
    
    bsl = dPEG_Ractive(i,1:tonset);
    dPEG_Ractive(i,:) = dPEG_Ractive(i,:) - nanmean(bsl);
    
    bsl = dPEG_Tactive(i,1:tonset);
    dPEG_Tactive(i,:) = dPEG_Tactive(i,:) - nanmean(bsl);
    
    bsl = dPEG_rCLT(i,1:tonset);
    dPEG_rCLT(i,:) = dPEG_rCLT(i,:) - nanmean(bsl);
    
    bsl = dPEG_tCLT(i,1:tonset);
    dPEG_tCLT(i,:) = dPEG_tCLT(i,:) - nanmean(bsl);
  end
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([dPEG_Rpassive(i,:) dPEG_Ractive(i,:) dPEG_Tpassive(i,:) dPEG_Tactive(i,:)]);
      fmaxP = fmax; fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([dPEG_Rpassive(i,:) dPEG_Tpassive(i,:)]);
      fmaxA = max([dPEG_Ractive(i,:) dPEG_Tactive(i,:)]);
    end
    
    dPEG_Rpassive(i,:) = dPEG_Rpassive(i,:) ./ fmaxP;
    dPEG_Tpassive(i,:) = dPEG_Tpassive(i,:) ./ fmaxP;
    dPEG_Ractive(i,:) = dPEG_Ractive(i,:)   ./ fmaxA;
    dPEG_Tactive(i,:) = dPEG_Tactive(i,:)   ./ fmaxA;
    
    dPEG_rCLT(i,:)     = dPEG_rCLT(i,:)       ./ fmaxP;
    dPEG_tCLT(i,:)     = dPEG_tCLT(i,:)       ./ fmaxP;
  end
  
  
end

% dPEG UP

for i=1:size(dPEG_RpassiveU, 1)
  % baseline mean removal
  
  tonset = 0.2 * dPEG_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = dPEG_RpassiveU(i,1:tonset);
    dPEG_RpassiveU(i,:) = dPEG_RpassiveU(i,:) - nanmean(bsl);
    
    bsl = dPEG_TpassiveU(i,1:tonset);
    dPEG_TpassiveU(i,:) = dPEG_TpassiveU(i,:) - nanmean(bsl);
    
    bsl = dPEG_RactiveU(i,1:tonset);
    dPEG_RactiveU(i,:) = dPEG_RactiveU(i,:) - nanmean(bsl);
    
    bsl = dPEG_TactiveU(i,1:tonset);
    dPEG_TactiveU(i,:) = dPEG_TactiveU(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([dPEG_RpassiveU(i,:) dPEG_RactiveU(i,:) dPEG_TpassiveU(i,:) dPEG_TactiveU(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([dPEG_RpassiveU(i,:) dPEG_TpassiveU(i,:)]);
      fmaxA = max([dPEG_RactiveU(i,:) dPEG_TactiveU(i,:)]);
    end
    
    dPEG_RpassiveU(i,:) = dPEG_RpassiveU(i,:) ./ fmaxP;
    dPEG_TpassiveU(i,:) = dPEG_TpassiveU(i,:) ./ fmaxP;
    dPEG_RactiveU(i,:) = dPEG_RactiveU(i,:)   ./ fmaxA;
    dPEG_TactiveU(i,:) = dPEG_TactiveU(i,:)   ./ fmaxA;
    
  end
    
end

% dPEG DOWN

for i=1:size(dPEG_RpassiveD, 1)
  % baseline mean removal
  
  tonset = 0.2 * dPEG_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = dPEG_RpassiveD(i,1:tonset);
    dPEG_RpassiveD(i,:) = dPEG_RpassiveD(i,:) - nanmean(bsl);
    
    bsl = dPEG_TpassiveD(i,1:tonset);
    dPEG_TpassiveD(i,:) = dPEG_TpassiveD(i,:) - nanmean(bsl);
    
    bsl = dPEG_RactiveD(i,1:tonset);
    dPEG_RactiveD(i,:) = dPEG_RactiveD(i,:) - nanmean(bsl);
    
    bsl = dPEG_TactiveD(i,1:tonset);
    dPEG_TactiveD(i,:) = dPEG_TactiveD(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([dPEG_RpassiveD(i,:) dPEG_RactiveD(i,:) dPEG_TpassiveD(i,:) dPEG_TactiveD(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([dPEG_RpassiveD(i,:) dPEG_TpassiveD(i,:)]);
      fmaxA = max([dPEG_RactiveD(i,:) dPEG_TactiveD(i,:)]);
    end
    
    dPEG_RpassiveD(i,:) = dPEG_RpassiveD(i,:) ./ fmaxP;
    dPEG_TpassiveD(i,:) = dPEG_TpassiveD(i,:) ./ fmaxP;
    dPEG_RactiveD(i,:) = dPEG_RactiveD(i,:)   ./ fmaxA;
    dPEG_TactiveD(i,:) = dPEG_TactiveD(i,:)   ./ fmaxA;
    
  end
    
end

if popnorm
  [dPEG_Rpassive, dPEG_Ractive, dPEG_Tpassive, dPEG_Tactive] = popmaxnorm(dPEG_Rpassive, dPEG_Ractive, dPEG_Tpassive, dPEG_Tactive);
  [dPEG_RpassiveU, dPEG_RactiveU, dPEG_TpassiveU, dPEG_TactiveU] = popmaxnorm(dPEG_RpassiveU, dPEG_RactiveU, dPEG_TpassiveU, dPEG_TactiveU);
  [dPEG_RpassiveD, dPEG_RactiveD, dPEG_TpassiveD, dPEG_TactiveD] = popmaxnorm(dPEG_RpassiveD, dPEG_RactiveD, dPEG_TpassiveD, dPEG_TactiveD);
end

if add_naive
  % naive (tango) data
  tonset = 0.2 * ndPEG_sr; % presilence was fixed at 0.2s
  
  for i=1:size(ndPEG_Rpassive, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = ndPEG_Rpassive(i,1:tonset);
      ndPEG_Rpassive(i,:) = ndPEG_Rpassive(i,:) - nanmean(bsl);
      
      bsl = ndPEG_Tpassive(i,1:tonset);
      ndPEG_Tpassive(i,:) = ndPEG_Tpassive(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([ndPEG_Rpassive(i,:) ndPEG_Tpassive(i,:)]);
      ndPEG_Rpassive(i,:) = ndPEG_Rpassive(i,:) ./ fmax;
      ndPEG_Tpassive(i,:) = ndPEG_Tpassive(i,:) ./ fmax;
    end
  end
  
  % click discrim going UP
  for i=1:size(ndPEG_RpassiveU, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = ndPEG_RpassiveU(i,1:tonset);
      ndPEG_RpassiveU(i,:) = ndPEG_RpassiveU(i,:) - nanmean(bsl);
      
      bsl = ndPEG_TpassiveU(i,1:tonset);
      ndPEG_TpassiveU(i,:) = ndPEG_TpassiveU(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([ndPEG_RpassiveU(i,:) ndPEG_TpassiveU(i,:)]);
      ndPEG_RpassiveU(i,:) = ndPEG_RpassiveU(i,:) ./ fmax;
      ndPEG_TpassiveU(i,:) = ndPEG_TpassiveU(i,:) ./ fmax;
    end
  end
  
  % click discrim going DOWN
  for i=1:size(ndPEG_RpassiveD, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = ndPEG_RpassiveD(i,1:tonset);
      ndPEG_RpassiveD(i,:) = ndPEG_RpassiveD(i,:) - nanmean(bsl);
      
      bsl = ndPEG_TpassiveD(i,1:tonset);
      ndPEG_TpassiveD(i,:) = ndPEG_TpassiveD(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([ndPEG_RpassiveD(i,:) ndPEG_TpassiveD(i,:)]);
      ndPEG_RpassiveD(i,:) = ndPEG_RpassiveD(i,:) ./ fmax;
      ndPEG_TpassiveD(i,:) = ndPEG_TpassiveD(i,:) ./ fmax;
    end
  end  
  
  if popnorm
    [ndPEG_Rpassive, ndPEG_Tpassive] = popmaxnorm(ndPEG_Rpassive, ndPEG_Tpassive);
    [ndPEG_RpassiveU, ndPEG_TpassiveU] = popmaxnorm(ndPEG_RpassiveU, ndPEG_TpassiveU);
    [ndPEG_RpassiveD, ndPEG_TpassiveD] = popmaxnorm(ndPEG_RpassiveD, ndPEG_TpassiveD);
  end
end


%% VPr

nbins = size(VPr_Rpassive, 2);

for i=1:size(VPr_Rpassive, 1)
  % baseline mean removal
  
  tonset = 0.2 * VPr_sr; % presilence was fixed at 0.2s
  %bsl = zeros(1,nbins).*NaN;
  
  VPrp_bsl = [VPrp_bsl nanmean([VPr_Rpassive(i,1:tonset) VPr_Tpassive(i,1:tonset)])];
  VPra_bsl = [VPra_bsl nanmean([VPr_Ractive(i,1:tonset)  VPr_Tactive(i,1:tonset)])];  
  
  if subtract_baseline
    bsl = VPr_Rpassive(i,1:tonset);
    VPr_Rpassive(i,:) = VPr_Rpassive(i,:) - nanmean(bsl);
    
    bsl = VPr_Tpassive(i,1:tonset);
    VPr_Tpassive(i,:) = VPr_Tpassive(i,:) - nanmean(bsl);
    
    bsl = VPr_Ractive(i,1:tonset);
    VPr_Ractive(i,:) = VPr_Ractive(i,:) - nanmean(bsl);
    
    bsl = VPr_Tactive(i,1:tonset);
    VPr_Tactive(i,:) = VPr_Tactive(i,:) - nanmean(bsl);
    
    bsl = VPr_rCLT(i,1:tonset);
    VPr_rCLT(i,:) = VPr_rCLT(i,:) - nanmean(bsl);
    
    bsl = VPr_tCLT(i,1:tonset);
    VPr_tCLT(i,:) = VPr_tCLT(i,:) - nanmean(bsl);
  end
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([VPr_Rpassive(i,:) VPr_Ractive(i,:) VPr_Tpassive(i,:) VPr_Tactive(i,:)]);
      fmaxP = fmax; fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([VPr_Rpassive(i,:) VPr_Tpassive(i,:)]);
      fmaxA = max([VPr_Ractive(i,:) VPr_Tactive(i,:)]);
    end
    
    VPr_Rpassive(i,:) = VPr_Rpassive(i,:) ./ fmaxP;
    VPr_Tpassive(i,:) = VPr_Tpassive(i,:) ./ fmaxP;
    VPr_Ractive(i,:) = VPr_Ractive(i,:)   ./ fmaxA;
    VPr_Tactive(i,:) = VPr_Tactive(i,:)   ./ fmaxA;
    
    VPr_rCLT(i,:)     = VPr_rCLT(i,:)       ./ fmaxP;
    VPr_tCLT(i,:)     = VPr_tCLT(i,:)       ./ fmaxP;
  end
  
  
end

% VPr UP

for i=1:size(VPr_RpassiveU, 1)
  % baseline mean removal
  
  tonset = 0.2 * VPr_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = VPr_RpassiveU(i,1:tonset);
    VPr_RpassiveU(i,:) = VPr_RpassiveU(i,:) - nanmean(bsl);
    
    bsl = VPr_TpassiveU(i,1:tonset);
    VPr_TpassiveU(i,:) = VPr_TpassiveU(i,:) - nanmean(bsl);
    
    bsl = VPr_RactiveU(i,1:tonset);
    VPr_RactiveU(i,:) = VPr_RactiveU(i,:) - nanmean(bsl);
    
    bsl = VPr_TactiveU(i,1:tonset);
    VPr_TactiveU(i,:) = VPr_TactiveU(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([VPr_RpassiveU(i,:) VPr_RactiveU(i,:) VPr_TpassiveU(i,:) VPr_TactiveU(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([VPr_RpassiveU(i,:) VPr_TpassiveU(i,:)]);
      fmaxA = max([VPr_RactiveU(i,:) VPr_TactiveU(i,:)]);
    end
    
    VPr_RpassiveU(i,:) = VPr_RpassiveU(i,:) ./ fmaxP;
    VPr_TpassiveU(i,:) = VPr_TpassiveU(i,:) ./ fmaxP;
    VPr_RactiveU(i,:) = VPr_RactiveU(i,:)   ./ fmaxA;
    VPr_TactiveU(i,:) = VPr_TactiveU(i,:)   ./ fmaxA;
    
  end
    
end

% VPr DOWN

for i=1:size(VPr_RpassiveD, 1)
  % baseline mean removal
  
  tonset = 0.2 * VPr_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = VPr_RpassiveD(i,1:tonset);
    VPr_RpassiveD(i,:) = VPr_RpassiveD(i,:) - nanmean(bsl);
    
    bsl = VPr_TpassiveD(i,1:tonset);
    VPr_TpassiveD(i,:) = VPr_TpassiveD(i,:) - nanmean(bsl);
    
    bsl = VPr_RactiveD(i,1:tonset);
    VPr_RactiveD(i,:) = VPr_RactiveD(i,:) - nanmean(bsl);
    
    bsl = VPr_TactiveD(i,1:tonset);
    VPr_TactiveD(i,:) = VPr_TactiveD(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([VPr_RpassiveD(i,:) VPr_RactiveD(i,:) VPr_TpassiveD(i,:) VPr_TactiveD(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([VPr_RpassiveD(i,:) VPr_TpassiveD(i,:)]);
      fmaxA = max([VPr_RactiveD(i,:) VPr_TactiveD(i,:)]);
    end
    
    VPr_RpassiveD(i,:) = VPr_RpassiveD(i,:) ./ fmaxP;
    VPr_TpassiveD(i,:) = VPr_TpassiveD(i,:) ./ fmaxP;
    VPr_RactiveD(i,:) = VPr_RactiveD(i,:)   ./ fmaxA;
    VPr_TactiveD(i,:) = VPr_TactiveD(i,:)   ./ fmaxA;
    
  end
    
end

if popnorm
  [VPr_Rpassive, VPr_Ractive, VPr_Tpassive, VPr_Tactive] = popmaxnorm(VPr_Rpassive, VPr_Ractive, VPr_Tpassive, VPr_Tactive);
  [VPr_RpassiveU, VPr_RactiveU, VPr_TpassiveU, VPr_TactiveU] = popmaxnorm(VPr_RpassiveU, VPr_RactiveU, VPr_TpassiveU, VPr_TactiveU);
  [VPr_RpassiveD, VPr_RactiveD, VPr_TpassiveD, VPr_TactiveD] = popmaxnorm(VPr_RpassiveD, VPr_RactiveD, VPr_TpassiveD, VPr_TactiveD);
end

if add_naive
  % naive (tango) data
  tonset = 0.2 * nVPr_sr; % presilence was fixed at 0.2s
  
  for i=1:size(nVPr_Rpassive, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = nVPr_Rpassive(i,1:tonset);
      nVPr_Rpassive(i,:) = nVPr_Rpassive(i,:) - nanmean(bsl);
      
      bsl = nVPr_Tpassive(i,1:tonset);
      nVPr_Tpassive(i,:) = nVPr_Tpassive(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([nVPr_Rpassive(i,:) nVPr_Tpassive(i,:)]);
      nVPr_Rpassive(i,:) = nVPr_Rpassive(i,:) ./ fmax;
      nVPr_Tpassive(i,:) = nVPr_Tpassive(i,:) ./ fmax;
    end
  end
  
  % click discrim going UP
  for i=1:size(nVPr_RpassiveU, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = nVPr_RpassiveU(i,1:tonset);
      nVPr_RpassiveU(i,:) = nVPr_RpassiveU(i,:) - nanmean(bsl);
      
      bsl = nVPr_TpassiveU(i,1:tonset);
      nVPr_TpassiveU(i,:) = nVPr_TpassiveU(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([nVPr_RpassiveU(i,:) nVPr_TpassiveU(i,:)]);
      nVPr_RpassiveU(i,:) = nVPr_RpassiveU(i,:) ./ fmax;
      nVPr_TpassiveU(i,:) = nVPr_TpassiveU(i,:) ./ fmax;
    end
  end
  
  % click discrim going DOWN
  for i=1:size(nVPr_RpassiveD, 1)
    % baseline mean removal
    if subtract_baseline
      bsl = nVPr_RpassiveD(i,1:tonset);
      nVPr_RpassiveD(i,:) = nVPr_RpassiveD(i,:) - nanmean(bsl);
      
      bsl = nVPr_TpassiveD(i,1:tonset);
      nVPr_TpassiveD(i,:) = nVPr_TpassiveD(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      fmax = max([nVPr_RpassiveD(i,:) nVPr_TpassiveD(i,:)]);
      nVPr_RpassiveD(i,:) = nVPr_RpassiveD(i,:) ./ fmax;
      nVPr_TpassiveD(i,:) = nVPr_TpassiveD(i,:) ./ fmax;
    end
  end  
  
  if popnorm
    [nVPr_Rpassive, nVPr_Tpassive] = popmaxnorm(nVPr_Rpassive, nVPr_Tpassive);
    [nVPr_RpassiveU, nVPr_TpassiveU] = popmaxnorm(nVPr_RpassiveU, nVPr_TpassiveU);
    [nVPr_RpassiveD, nVPr_TpassiveD] = popmaxnorm(nVPr_RpassiveD, nVPr_TpassiveD);
  end
end

%% PFC

nbins = size(PFC_Rpassive, 2);

for i=1:size(PFC_Rpassive, 1)
  % baseline mean removal
  
  tonset = 0.2 * PFC_sr; % presilence was fixed at 0.2s
  %bsl = zeros(1,nbins).*NaN;
  
  PFCp_bsl = [PFCp_bsl nanmean([PFC_Rpassive(i,1:tonset) PFC_Tpassive(i,1:tonset)])];
  PFCa_bsl = [PFCa_bsl nanmean([PFC_Ractive(i,1:tonset)  PFC_Tactive(i,1:tonset)])];  
  
  if subtract_baseline
    bsl = PFC_Rpassive(i,1:tonset);
    PFC_Rpassive(i,:) = PFC_Rpassive(i,:) - nanmean(bsl);
    
    bsl = PFC_Tpassive(i,1:tonset);
    PFC_Tpassive(i,:) = PFC_Tpassive(i,:) - nanmean(bsl);
    
    bsl = PFC_Ractive(i,1:tonset);
    PFC_Ractive(i,:) = PFC_Ractive(i,:) - nanmean(bsl);
    
    bsl = PFC_Tactive(i,1:tonset);
    PFC_Tactive(i,:) = PFC_Tactive(i,:) - nanmean(bsl);
    
    bsl = PFC_rCLT(i,1:tonset);
    PFC_rCLT(i,:) = PFC_rCLT(i,:) - nanmean(bsl);
    
    bsl = PFC_tCLT(i,1:tonset);
    PFC_tCLT(i,:) = PFC_tCLT(i,:) - nanmean(bsl);
  end
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([PFC_Rpassive(i,:) PFC_Ractive(i,:) PFC_Tpassive(i,:) PFC_Tactive(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([PFC_Rpassive(i,:) PFC_Tpassive(i,:)]);
      fmaxA = max([PFC_Ractive(i,:) PFC_Tactive(i,:)]);
    end
    
    PFC_Rpassive(i,:) = PFC_Rpassive(i,:) ./ fmax;
    PFC_Tpassive(i,:) = PFC_Tpassive(i,:) ./ fmax;
    PFC_Ractive(i,:) = PFC_Ractive(i,:)   ./ fmax;
    PFC_Tactive(i,:) = PFC_Tactive(i,:)   ./ fmax;
    
    PFC_rCLT(i,:)     = PFC_rCLT(i,:)       ./ fmax;
    PFC_tCLT(i,:)     = PFC_tCLT(i,:)       ./ fmax;
  end
  
end

% PFC UP

for i=1:size(PFC_RpassiveU, 1)
  % baseline mean removal
  
  tonset = 0.2 * PFC_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = PFC_RpassiveU(i,1:tonset);
    PFC_RpassiveU(i,:) = PFC_RpassiveU(i,:) - nanmean(bsl);
    
    bsl = PFC_TpassiveU(i,1:tonset);
    PFC_TpassiveU(i,:) = PFC_TpassiveU(i,:) - nanmean(bsl);
    
    bsl = PFC_RactiveU(i,1:tonset);
    PFC_RactiveU(i,:) = PFC_RactiveU(i,:) - nanmean(bsl);
    
    bsl = PFC_TactiveU(i,1:tonset);
    PFC_TactiveU(i,:) = PFC_TactiveU(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([PFC_RpassiveU(i,:) PFC_RactiveU(i,:) PFC_TpassiveU(i,:) PFC_TactiveU(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([PFC_RpassiveU(i,:) PFC_TpassiveU(i,:)]);
      fmaxA = max([PFC_RactiveU(i,:) PFC_TactiveU(i,:)]);
    end
    
    PFC_RpassiveU(i,:) = PFC_RpassiveU(i,:) ./ fmaxP;
    PFC_TpassiveU(i,:) = PFC_TpassiveU(i,:) ./ fmaxP;
    PFC_RactiveU(i,:) = PFC_RactiveU(i,:)   ./ fmaxA;
    PFC_TactiveU(i,:) = PFC_TactiveU(i,:)   ./ fmaxA;
    
  end
    
end

% PFC DOWN

for i=1:size(PFC_RpassiveD, 1)
  % baseline mean removal
  
  tonset = 0.2 * PFC_sr; % presilence was fixed at 0.2s
  
  if subtract_baseline
    bsl = PFC_RpassiveD(i,1:tonset);
    PFC_RpassiveD(i,:) = PFC_RpassiveD(i,:) - nanmean(bsl);
    
    bsl = PFC_TpassiveD(i,1:tonset);
    PFC_TpassiveD(i,:) = PFC_TpassiveD(i,:) - nanmean(bsl);
    
    bsl = PFC_RactiveD(i,1:tonset);
    PFC_RactiveD(i,:) = PFC_RactiveD(i,:) - nanmean(bsl);
    
    bsl = PFC_TactiveD(i,1:tonset);
    PFC_TactiveD(i,:) = PFC_TactiveD(i,:) - nanmean(bsl);
    
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 0
      fmax = max([PFC_RpassiveD(i,:) PFC_RactiveD(i,:) PFC_TpassiveD(i,:) PFC_TactiveD(i,:)]);
      fmaxP = fmax;
      fmaxA = fmax;
    elseif indep_norm == 1
      fmaxP = max([PFC_RpassiveD(i,:) PFC_TpassiveD(i,:)]);
      fmaxA = max([PFC_RactiveD(i,:) PFC_TactiveD(i,:)]);
    end
    
    PFC_RpassiveD(i,:) = PFC_RpassiveD(i,:) ./ fmaxP;
    PFC_TpassiveD(i,:) = PFC_TpassiveD(i,:) ./ fmaxP;
    PFC_RactiveD(i,:) = PFC_RactiveD(i,:)   ./ fmaxA;
    PFC_TactiveD(i,:) = PFC_TactiveD(i,:)   ./ fmaxA;
    
  end
    
end

if popnorm
  [PFC_Rpassive, PFC_Ractive, PFC_Tpassive, PFC_Tactive] = popmaxnorm(PFC_Rpassive, PFC_Ractive, PFC_Tpassive, PFC_Tactive);
  [PFC_RpassiveU, PFC_RactiveU, PFC_TpassiveU, PFC_TactiveU] = popmaxnorm(PFC_RpassiveU, PFC_RactiveU, PFC_TpassiveU, PFC_TactiveU);
  [PFC_RpassiveD, PFC_RactiveD, PFC_TpassiveD, PFC_TactiveD] = popmaxnorm(PFC_RpassiveD, PFC_RactiveD, PFC_TpassiveD, PFC_TactiveD);
end


%% svd_dPEG: subset of dPEG neurons that were used in the neuron paper (2014)
if svd_subset
  nbins = size(svd_dPEG_Rpassive, 2);
  
  for i=1:size(svd_dPEG_Rpassive, 1)
    % baseline mean removal
    
    tonset = 0.2 * svd_dPEG_sr; % presilence was fixed at 0.2s
    %bsl = zeros(1,nbins).*NaN;
    
    
    if subtract_baseline
      bsl = svd_dPEG_Rpassive(i,1:tonset);
      svd_dPEG_Rpassive(i,:) = svd_dPEG_Rpassive(i,:) - mean(bsl);
      
      bsl = svd_dPEG_Tpassive(i,1:tonset);
      svd_dPEG_Tpassive(i,:) = svd_dPEG_Tpassive(i,:) - nanmean(bsl);
      
      bsl = svd_dPEG_Ractive(i,1:tonset);
      svd_dPEG_Ractive(i,:) = svd_dPEG_Ractive(i,:) - nanmean(bsl);
      
      bsl = svd_dPEG_Tactive(i,1:tonset);
      svd_dPEG_Tactive(i,:) = svd_dPEG_Tactive(i,:) - nanmean(bsl);
      
      bsl = svd_dPEG_rCLT(i,1:tonset);
      svd_dPEG_rCLT(i,:) = svd_dPEG_rCLT(i,:) - nanmean(bsl);
      
      bsl = svd_dPEG_tCLT(i,1:tonset);
      svd_dPEG_tCLT(i,:) = svd_dPEG_tCLT(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      if indep_norm == 0
        fmax = max([svd_dPEG_Rpassive(i,:) svd_dPEG_Ractive(i,:) svd_dPEG_Tpassive(i,:) svd_dPEG_Tactive(i,:)]);
        fmaxP = fmax;
        fmaxA = fmax;
      elseif indep_norm == 1
        fmaxP = max([svd_dPEG_Rpassive(i,:) svd_dPEG_Tpassive(i,:)]);
        fmaxA = max([svd_dPEG_Ractive(i,:) svd_dPEG_Tactive(i,:)]);
      end
      
      svd_dPEG_Rpassive(i,:) = svd_dPEG_Rpassive(i,:) ./ fmax;
      svd_dPEG_Tpassive(i,:) = svd_dPEG_Tpassive(i,:) ./ fmax;
      svd_dPEG_Ractive(i,:) = svd_dPEG_Ractive(i,:)   ./ fmax;
      svd_dPEG_Tactive(i,:) = svd_dPEG_Tactive(i,:)   ./ fmax;
      
      svd_dPEG_rCLT(i,:)     = svd_dPEG_rCLT(i,:)       ./ fmax;
      svd_dPEG_tCLT(i,:)     = svd_dPEG_tCLT(i,:)       ./ fmax;
      
    end
    
    
  end
  
  if popnorm
    [svd_dPEG_Rpassive, svd_dPEG_Ractive, svd_dPEG_Tpassive, svd_dPEG_Tactive] = popmaxnorm(svd_dPEG_Rpassive, svd_dPEG_Ractive, svd_dPEG_Tpassive, svd_dPEG_Tactive);
  end
  
  svd_dPEG_t(end) = []; % there was an extra bin at the end, check why
end

%% Selection rules

%%---A1
% selection of cells based on target positive *onset* responses (tpr) during behavior
clkon = zeros(length(A1_params),1);
for i = 1: length(A1_params)
  clkon(i) = ceil((A1_params{i}.TrialObject.ReferenceHandle.TorcDuration + 0.4) .* A1_sr);
end


if only_pos_tar == 1
  if pos_tar_act == 1
    A1tpr = [find(mean(A1_Tactive(:,scatwin),2)>0)];
  elseif pos_tar_act == 0
    A1tpr = [find(mean(A1_Tpassive(:,scatwin),2)>0)];
  end
  
elseif only_pos_tar ==-1
  if pos_tar_act == 1
    A1tpr = [find(mean(A1_Tactive(:,scatwin),2)<=0)];
  elseif pos_tar_act == 0
    A1tpr = [find(mean(A1_Tpassive(:,scatwin),2)<=0)];
  end
  
elseif only_pos_tar == 0
  A1tpr = 1:size(A1_Tactive, 1);
end

% leave only 1 target files
if limit_single_tar == 1
  aa = zeros(length(A1tpr),1) .* NaN;
  for i = 1:length(A1tpr)
    j=A1tpr(i);
    if length(A1_params{j}.TrialObject.TargetHandle.Frequencies) == 1
      aa(i) = j;
    end
  end
  aa(isnan(aa))=[];
  A1tpr = aa;
  
elseif limit_single_tar == -1
  aa = zeros(length(A1tpr),1) .* NaN;
  for i = 1:length(A1tpr)
    j=A1tpr(i);
    if length(A1_params{j}.TrialObject.TargetHandle.Frequencies) > 1
      aa(i) = j;
    end
  end
  aa(isnan(aa))=[];
  A1tpr = aa;
end


%%---dPEG
if svd_subset == 0
  
  if only_pos_tar == 1
    if pos_tar_act == 0
      dPEGtpr = [find(mean(dPEG_Tpassive(:,50:60),2)>0)];
    elseif pos_tar_act == 1
      dPEGtpr = [find(mean(dPEG_Tactive(:,scatwin),2)>0)];
    end
    %dPEGtpr = find(mean([dPEG_Tactive(:,50:60) dPEG_Tpassive(:,50:60)],2)>0);
  elseif only_pos_tar == -1
    if pos_tar_act == 1
      dPEGtpr = [find(mean(dPEG_Tactive(:,scatwin),2)<=0)];
    elseif pos_tar_act == 0
      dPEGtpr = [find(mean(dPEG_Tpassive(:,scatwin),2)<=0)];
    end
    
  elseif only_pos_tar == 0
    dPEGtpr = 1:size(dPEG_Tactive, 1);
  end
  
  
  % leave only 1 target files
  if limit_single_tar == 1
    aa = zeros(length(dPEGtpr),1) .* NaN;
    for i = 1:length(dPEGtpr)
      j=dPEGtpr(i);
      if length(dPEG_params{j}.TrialObject.TargetHandle.Frequencies) == 1
        aa(i) = j;
      end
    end
    aa(isnan(aa))=[];
    dPEGtpr = aa;
    
  elseif limit_single_tar == -1
    aa = zeros(length(dPEGtpr),1) .* NaN;
    for i = 1:length(dPEGtpr)
      j=dPEGtpr(i);
      if length(dPEG_params{j}.TrialObject.TargetHandle.Frequencies) > 1
        aa(i) = j;
      end
    end
    aa(isnan(aa))=[];
    dPEGtpr = aa;
    
  end
end


%%--svd_dPEG
if svd_subset == 1
  if only_pos_tar == 1
    if pos_tar_act == 0
      dPEGtpr = [find(mean(svd_dPEG_Tpassive(:,50:60),2)>0)];
    elseif pos_tar_act == 1
      dPEGtpr = [find(mean(svd_dPEG_Tactive(:,50:60),2)>0)];
    end
    %dPEGtpr = find(mean([svd_dPEG_Tactive(:,50:60) svd_dPEG_Tpassive(:,50:60)],2)>0);
  elseif only_pos_tar == -1
    if pos_tar_act == 1
      dPEGtpr = [find(mean(svd_dPEG_Tactive(:,scatwin),2)<=0)];
    elseif pos_tar_act == 0
      dPEGtpr = [find(mean(svd_dPEG_Tpassive(:,scatwin),2)<=0)];
    end
    
  elseif only_pos_tar == 0
    dPEGtpr = 1:size(svd_dPEG_Tactive, 1);
  end
  
  % leave only 1 target files
  if limit_single_tar == 1
    aa = zeros(length(dPEGtpr),1) .* NaN;
    for i = 1:length(dPEGtpr)
      j=dPEGtpr(i);
      if length(svd_dPEG_params{j}.TrialObject.TargetHandle.Frequencies) == 1
        aa(i) = j;
      end
    end
    aa(isnan(aa))=[];
    dPEGtpr = aa;
    
  elseif limit_single_tar == -1
    aa = zeros(length(dPEGtpr),1) .* NaN;
    for i = 1:length(dPEGtpr)
      j=dPEGtpr(i);
      if length(svd_dPEG_params{j}.TrialObject.TargetHandle.Frequencies) > 1
        aa(i) = j;
      end
    end
    aa(isnan(aa))=[];
    dPEGtpr = aa;
  end
end

%%--VPr
if only_pos_tar ==1
  if pos_tar_act == 1
    VPrtpr = [find(mean(VPr_Tactive(:,scatwin),2)>0)];
  elseif pos_tar_act == 0
    VPrtpr = [find(mean(VPr_Tpassive(:,scatwin),2)>0)];
  end
  
elseif only_pos_tar == -1
  if pos_tar_act == 1
    VPrtpr = [find(mean(VPr_Tactive(:,scatwin),2)<=0)];
  elseif pos_tar_act == 0
    VPrtpr = [find(mean(VPr_Tpassive(:,scatwin),2)<=0)];
  end
  
elseif only_pos_tar == 0
  VPrtpr = 1:size(VPr_Tactive,1);
end

% leave only 1 target files
if limit_single_tar == 1
  aa = zeros(length(VPrtpr),1) .* NaN;
  for i = 1:length(VPrtpr)
    j=VPrtpr(i);
    if length(VPr_params{j}.TrialObject.TargetHandle.Frequencies) == 1
      aa(i) = j;
    end
  end
  aa(isnan(aa))=[];
  VPrtpr = aa;
  
elseif limit_single_tar == -1
  aa = zeros(length(VPrtpr),1) .* NaN;
  for i = 1:length(VPrtpr)
    j=VPrtpr(i);
    if length(VPr_params{j}.TrialObject.TargetHandle.Frequencies) > 1
      aa(i) = j;
    end
  end
  aa(isnan(aa))=[];
  VPrtpr = aa;
end


%%--PFC
if only_pos_tar == 1
  if pos_tar_act == 1
    PFCtpr = [find(mean(PFC_Tactive(:,scatwin),2)>0)];
  elseif pos_tar_act == 0
    PFCtpr = [find(mean(PFC_Tpassive(:,scatwin),2)>0)];
  end
  
elseif only_pos_tar == -1
  if pos_tar_act == 1
    PFCtpr = [find(mean(PFC_Tactive(:,scatwin),2)<=0)];
  elseif pos_tar_act == 0
    PFCtpr = [find(mean(PFC_Tpassive(:,scatwin),2)<=0)];
  end
  
elseif only_pos_tar == 0
  PFCtpr = 1:size(PFC_Tactive, 1);
end

% leave only 1 target files
if limit_single_tar == 1
  aa = zeros(length(PFCtpr),1) .* NaN;
  for i = 1:length(PFCtpr)
    j=PFCtpr(i);
    if length(PFC_params{j}.TrialObject.TargetHandle.Frequencies) == 1
      aa(i) = j;
    end
  end
  aa(isnan(aa))=[];
  PFCtpr = aa;
  
elseif limit_single_tar == -1
  aa = zeros(length(PFCtpr),1) .* NaN;
  for i = 1:length(PFCtpr)
    j=PFCtpr(i);
    if length(PFC_params{j}.TrialObject.TargetHandle.Frequencies) > 1
      aa(i) = j;
    end
  end
  aa(isnan(aa))=[];
  PFCtpr = aa;
end

% simple fix to time axes
A1_t = A1_t - 0.4;
dPEG_t = dPEG_t - 0.4;
VPr_t = VPr_t - 0.4;
PFC_t = PFC_t - 0.4;

%% removal of last bins for cleaning the plot...
A1_Rpassive(:,97:end) = NaN;
A1_Tpassive(:,97:end) = NaN;
A1_Ractive(:,97:end) = NaN;
A1_Tactive(:,97:end) = NaN;


%% mean + SEM PSTHs

[mA1rp seA1rp] = jackmeanerr(A1_Rpassive(A1tpr, :), 20);
[mA1ra seA1ra] = jackmeanerr(A1_Ractive(A1tpr, :), 20);
[mA1tp seA1tp] = jackmeanerr(A1_Tpassive(A1tpr, :), 20);
[mA1ta seA1ta] = jackmeanerr(A1_Tactive(A1tpr, :), 20);

if svd_subset == 0
  [mdPEGrp sedPEGrp] = jackmeanerr(dPEG_Rpassive(dPEGtpr, :), 20);
  [mdPEGra sedPEGra] = jackmeanerr(dPEG_Ractive (dPEGtpr, :), 20);
  [mdPEGtp sedPEGtp] = jackmeanerr(dPEG_Tpassive(dPEGtpr, :), 20);
  [mdPEGta sedPEGta] = jackmeanerr(dPEG_Tactive (dPEGtpr, :), 20);
else
  [mdPEGrp sedPEGrp] = jackmeanerr(svd_dPEG_Rpassive(dPEGtpr, :), 20);
  [mdPEGra sedPEGra] = jackmeanerr(svd_dPEG_Ractive (dPEGtpr, :), 20);
  [mdPEGtp sedPEGtp] = jackmeanerr(svd_dPEG_Tpassive(dPEGtpr, :), 20);
  [mdPEGta sedPEGta] = jackmeanerr(svd_dPEG_Tactive (dPEGtpr, :), 20);
  dPEG_t = svd_dPEG_t;
end

[mVPrrp seVPrrp] = jackmeanerr(VPr_Rpassive(VPrtpr, :), 20);
[mVPrra seVPrra] = jackmeanerr(VPr_Ractive (VPrtpr, :), 20);
[mVPrtp seVPrtp] = jackmeanerr(VPr_Tpassive(VPrtpr, :), 20);
[mVPrta seVPrta] = jackmeanerr(VPr_Tactive (VPrtpr, :), 20);

[mPFCrp sePFCrp] = jackmeanerr(PFC_Rpassive(PFCtpr, :), 20);
[mPFCra sePFCra] = jackmeanerr(PFC_Ractive (PFCtpr, :), 20);
[mPFCtp sePFCtp] = jackmeanerr(PFC_Tpassive(PFCtpr, :), 20);
[mPFCta sePFCta] = jackmeanerr(PFC_Tactive (PFCtpr, :), 20);

% DISCRIMINATION FROM LOW TO HIGH (GOING UP)
[mA1rpU seA1rpU] = jackmeanerr(A1_RpassiveU, 20);
[mA1raU seA1raU] = jackmeanerr(A1_RactiveU, 20);
[mA1tpU seA1tpU] = jackmeanerr(A1_TpassiveU, 20);
[mA1taU seA1taU] = jackmeanerr(A1_TactiveU, 20);

if svd_subset == 0
  [mdPEGrpU sedPEGrpU] = jackmeanerr(dPEG_RpassiveU, 20);
  [mdPEGraU sedPEGraU] = jackmeanerr(dPEG_RactiveU , 20);
  [mdPEGtpU sedPEGtpU] = jackmeanerr(dPEG_TpassiveU, 20);
  [mdPEGtaU sedPEGtaU] = jackmeanerr(dPEG_TactiveU , 20);
else
  [mdPEGrpU sedPEGrpU] = jackmeanerr(svd_dPEG_RpassiveU, 20);
  [mdPEGraU sedPEGraU] = jackmeanerr(svd_dPEG_RactiveU , 20);
  [mdPEGtpU sedPEGtpU] = jackmeanerr(svd_dPEG_TpassiveU, 20);
  [mdPEGtaU sedPEGtaU] = jackmeanerr(svd_dPEG_TactiveU , 20);
end

[mVPrrpU seVPrrpU] = jackmeanerr(VPr_RpassiveU, 20);
[mVPrraU seVPrraU] = jackmeanerr(VPr_RactiveU , 20);
[mVPrtpU seVPrtpU] = jackmeanerr(VPr_TpassiveU, 20);
[mVPrtaU seVPrtaU] = jackmeanerr(VPr_TactiveU , 20);

[mPFCrpU sePFCrpU] = jackmeanerr(PFC_RpassiveU, 20);
[mPFCraU sePFCraU] = jackmeanerr(PFC_RactiveU , 20);
[mPFCtpU sePFCtpU] = jackmeanerr(PFC_TpassiveU, 20);
[mPFCtaU sePFCtaU] = jackmeanerr(PFC_TactiveU , 20);

% DISCRIMINATION FROM HIGH TO LOW (GOING DOWN)
[mA1rpD seA1rpD] = jackmeanerr(A1_RpassiveD, 20);
[mA1raD seA1raD] = jackmeanerr(A1_RactiveD, 20);
[mA1tpD seA1tpD] = jackmeanerr(A1_TpassiveD, 20);
[mA1taD seA1taD] = jackmeanerr(A1_TactiveD, 20);

if svd_subset == 0
  [mdPEGrpD sedPEGrpD] = jackmeanerr(dPEG_RpassiveD, 20);
  [mdPEGraD sedPEGraD] = jackmeanerr(dPEG_RactiveD , 20);
  [mdPEGtpD sedPEGtpD] = jackmeanerr(dPEG_TpassiveD, 20);
  [mdPEGtaD sedPEGtaD] = jackmeanerr(dPEG_TactiveD , 20);
else
  [mdPEGrpD sedPEGrpD] = jackmeanerr(svd_dPEG_RpassiveD, 20);
  [mdPEGraD sedPEGraD] = jackmeanerr(svd_dPEG_RactiveD , 20);
  [mdPEGtpD sedPEGtpD] = jackmeanerr(svd_dPEG_TpassiveD, 20);
  [mdPEGtaD sedPEGtaD] = jackmeanerr(svd_dPEG_TactiveD , 20);
end

[mVPrrpD seVPrrpD] = jackmeanerr(VPr_RpassiveD, 20);
[mVPrraD seVPrraD] = jackmeanerr(VPr_RactiveD , 20);
[mVPrtpD seVPrtpD] = jackmeanerr(VPr_TpassiveD, 20);
[mVPrtaD seVPrtaD] = jackmeanerr(VPr_TactiveD , 20);

[mPFCrpD sePFCrpD] = jackmeanerr(PFC_RpassiveD, 20);
[mPFCraD sePFCraD] = jackmeanerr(PFC_RactiveD , 20);
[mPFCtpD sePFCtpD] = jackmeanerr(PFC_TpassiveD, 20);
[mPFCtaD sePFCtaD] = jackmeanerr(PFC_TactiveD , 20);


% naive
[mnA1rp senA1rp] = jackmeanerr(nA1_Rpassive, 20);
[mnA1tp senA1tp] = jackmeanerr(nA1_Tpassive, 20);
[mndPEGrp sendPEGrp] = jackmeanerr(ndPEG_Rpassive, 20);
[mndPEGtp sendPEGtp] = jackmeanerr(ndPEG_Tpassive, 20);
[mnVPrrp senVPrrp] = jackmeanerr(nVPr_Rpassive, 20);
[mnVPrtp senVPrtp] = jackmeanerr(nVPr_Tpassive, 20);

[mnA1rpU senA1rpU] = jackmeanerr(nA1_RpassiveU, 20);
[mnA1tpU senA1tpU] = jackmeanerr(nA1_TpassiveU, 20);
[mndPEGrpU sendPEGrpU] = jackmeanerr(ndPEG_RpassiveU, 20);
[mndPEGtpU sendPEGtpU] = jackmeanerr(ndPEG_TpassiveU, 20);
[mnVPrrpU senVPrrpU] = jackmeanerr(nVPr_RpassiveU, 20);
[mnVPrtpU senVPrtpU] = jackmeanerr(nVPr_TpassiveU, 20);

[mnA1rpD senA1rpD] = jackmeanerr(nA1_RpassiveD, 20);
[mnA1tpD senA1tpD] = jackmeanerr(nA1_TpassiveD, 20);
[mndPEGrpD sendPEGrpD] = jackmeanerr(ndPEG_RpassiveD, 20);
[mndPEGtpD sendPEGtpD] = jackmeanerr(ndPEG_TpassiveD, 20);
[mnVPrrpD senVPrrpD] = jackmeanerr(nVPr_RpassiveD, 20);
[mnVPrtpD senVPrtpD] = jackmeanerr(nVPr_TpassiveD, 20);





%% plot 1: passive vs. active

figure;
set(gcf, 'Position', [200 52 509 739]);
if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.4 1.3];
elseif popnorm == 0 && normalize_each == 0
  yrange = [-5 15];
  
end

% A1

subplot(4,2,1);
shadedErrorBar(A1_t, mA1rp, seA1rp, '--b', 0);
hold on;
shadedErrorBar(A1_t, mA1ra, seA1ra, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1: Reference pass vs act');

subplot(4,2,2);
shadedErrorBar(A1_t, mA1tp, seA1tp, '--r', 0);
hold on;
shadedErrorBar(A1_t, mA1ta, seA1ta, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(A1_cells_used)))]);

% dPEG

subplot(4,2,3);
shadedErrorBar(dPEG_t, mdPEGrp, sedPEGrp, '--b', 0);
hold on;
shadedErrorBar(dPEG_t, mdPEGra, sedPEGra, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG: Reference pass vs act');

subplot(4,2,4);
shadedErrorBar(dPEG_t, mdPEGtp, sedPEGtp, '--r', 0);
hold on;
shadedErrorBar(dPEG_t, mdPEGta, sedPEGta, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(dPEG_cells_used)))]);




% VPr
subplot(4,2,5);
shadedErrorBar(VPr_t, mVPrrp, seVPrrp, '--b', 0);
hold on;
shadedErrorBar(VPr_t, mVPrra, seVPrra, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr: Reference pass vs act');

subplot(4,2,6);
shadedErrorBar(VPr_t, mVPrtp, seVPrtp, '--r', 0);
hold on;
shadedErrorBar(VPr_t, mVPrta, seVPrta, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(VPr_cells_used)))]);

% dlFC
subplot(4,2,7);
shadedErrorBar(PFC_t, mPFCrp, sePFCrp, '--b', 0);
hold on;
shadedErrorBar(PFC_t, mPFCra, sePFCra, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(ylim);
title('dlFC: Reference pass vs act');

subplot(4,2,8);
shadedErrorBar(PFC_t, mPFCtp, sePFCtp, '--r', 0);
hold on;
shadedErrorBar(PFC_t, mPFCta, sePFCta, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dlFC: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(PFC_cells_used)))]);

suptitle('CRD, passive vs. active');

set(gcf, 'PaperPositionMode', 'auto');
print('01', '-dpdf');

%% plot 2: Reference vs Target

figure;
set(gcf, 'Position', [300 52 509 739]);
set(gcf, 'Color', 'w');
if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.4 1.3];
elseif ~popnorm && ~normalize_each
  yrange = [-20 20];
end

% A1
subplot(4,2,1);
hp1 = shadedErrorBar(A1_t, mA1rp, seA1rp, '--b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1=shadedErrorBar(A1_t, mA1tp, seA1tp, '--r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('\bfPassive', 'FontSize', 14);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');

subplot(4,2,2);
hp1 = shadedErrorBar(A1_t, mA1ra, seA1ra, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(A1_t, mA1ta, seA1ta, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('\bfBehavior', 'FontSize', 14);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(A1tpr))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);

% dPEG
if svd_subset == 0
  subplot(4,2,3);
  hp1 = shadedErrorBar(dPEG_t, mdPEGrp, sedPEGrp, '--b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(dPEG_t, mdPEGtp, sedPEGtp, '--r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  
  subplot(4,2,4);
  hp1 = shadedErrorBar(dPEG_t, mdPEGra, sedPEGra, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(dPEG_t, mdPEGta, sedPEGta, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(dPEGtpr))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);
end

%svd dPEG
if svd_subset == 1
  subplot(4,2,3);
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGrp, sedPEGrp, '--b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, '--r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  
  subplot(4,2,4);
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGra, sedPEGra, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGta, sedPEGta, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(dPEGtpr))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);
end

% VPr
subplot(4,2,5);
hp1 = shadedErrorBar(VPr_t, mVPrrp, seVPrrp, '--b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(VPr_t, mVPrtp, seVPrtp, '--r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');

subplot(4,2,6);
hp1 = shadedErrorBar(VPr_t, mVPrra, seVPrra, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(VPr_t, mVPrta, seVPrta, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(VPrtpr))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);

% dlFC
subplot(4,2,7);
rh = shadedErrorBar(PFC_t, mPFCrp, sePFCrp, '--b', 0);
set(rh.edge, 'LineStyle', 'none');
set(rh.mainLine, 'LineWidth', 2);
hold on;
th = shadedErrorBar(PFC_t, mPFCtp, sePFCtp, '--r',   0);
set(th.edge, 'LineStyle', 'none');
set(th.mainLine, 'LineWidth', 2);
legend([rh.mainLine th.mainLine], 'Reference', 'Target');
legend BOXOFF;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');
xlabel('Time (s)');

subplot(4,2,8);
hp1 = shadedErrorBar(PFC_t, mPFCra, sePFCra, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(PFC_t, mPFCta, sePFCta, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(PFCtpr))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);
xlabel('Time (s)');

%suptitle('Target detection, reference vs. target');

set(gcf, 'PaperPositionMode', 'auto');
print('02', '-dpdf');



if latstat
  keyboard;
  alfa = 0.05;
  twin = 73; % do test from bin t0 to 73 (during click train only)
  t0 = 52; % 52 is t=1.3s bin (click train onset)
  
  A1p_tar   = zeros(twin-t0, 1);
  dPEGp_tar = zeros(twin-t0, 1);
  VPrp_tar  = zeros(twin-t0, 1);
  PFCp_tar  = zeros(twin-t0, 1);
  
  A1p_ref   = zeros(twin-t0, 1);
  dPEGp_ref = zeros(twin-t0, 1);
  VPrp_ref  = zeros(twin-t0, 1);
  PFCp_ref  = zeros(twin-t0, 1);
  
  % target active vs target passive
  cnt = 1;
  for n = t0 : twin
    A1p_tar(cnt)   = signrank(A1_Tpassive(:,n), A1_Tactive(:,n), 'alpha', alfa);
    dPEGp_tar(cnt) = signrank(dPEG_Tpassive(:,n), dPEG_Tactive(:,n), 'alpha', alfa);
    VPrp_tar(cnt)  = signrank(VPr_Tpassive(:,n), VPr_Tactive(:,n), 'alpha', alfa);
    PFCp_tar(cnt)  = signrank(PFC_Tpassive(:,n), PFC_Tactive(:,n), 'alpha', alfa);
    cnt = cnt +1;
  end
  
    % reference active vs reference passive
    cnt = 1;
  for n = t0 : twin
    A1p_ref(cnt)   = signrank(A1_Rpassive(:,n), A1_Ractive(:,n), 'alpha', alfa);
    dPEGp_ref(cnt) = signrank(dPEG_Rpassive(:,n), dPEG_Ractive(:,n), 'alpha', alfa);
    VPrp_ref(cnt)  = signrank(VPr_Rpassive(:,n), VPr_Ractive(:,n), 'alpha', alfa);
    PFCp_ref(cnt)  = signrank(PFC_Rpassive(:,n), PFC_Ractive(:,n), 'alpha', alfa);
    cnt = cnt +1;
  end
  
  A1sig_ref   = findstr(A1p_ref(t0:end)', [1,1,1]);
  dPEGsig_ref = findstr(dPEGp_ref(t0:end)', [1,1,1]);
  VPrsig_ref  = findstr(VPrp_ref(t0:end)', [1,1,1]);
  PFCsig_ref  = findstr(PFCp_ref(t0:end)', [1,1,1]);
  
  A1sig_tar   = findstr(A1p_tar(t0:end)', [1,1,1]);
  dPEGsig_tar = findstr(dPEGp_tar(t0:end)', [1,1,1]);
  VPrsig_tar  = findstr(VPrp_tar(t0:end)', [1,1,1]);
  PFCsig_tar  = findstr(PFCp_tar(t0:end)', [1,1,1]);


[A1p_ref<(alfa/twin), dPEGp_ref<(alfa/twin), VPrp_ref<(alfa/twin), PFCp_ref< (alfa/twin)]
[A1p_tar<(alfa/(twin-t0)), dPEGp_tar<(alfa/(twin-t0)), VPrp_tar<(alfa/(twin-t0)), PFCp_tar< (alfa/(twin-t0))]


[A1p_ref<alfa, dPEGp_ref<alfa, VPrp_ref<alfa, PFCp_ref< alfa]
[A1p_tar<alfa, dPEGp_tar<alfa, VPrp_tar<alfa, PFCp_tar< alfa]

%no correction
[A1p_ref<alfa | A1p_tar<alfa, dPEGp_ref<alfa | dPEGp_tar<alfa, VPrp_ref<alfa | VPrp_tar<alfa, PFCp_ref< alfa | PFCp_tar< alfa]

%bonferroni
[A1p_ref<(alfa/(twin-t0)) | A1p_tar<(alfa/(twin-t0)), dPEGp_ref<(alfa/(twin-t0)) | dPEGp_tar<(alfa/(twin-t0)), VPrp_ref<(alfa/(twin-t0)) | VPrp_tar<(alfa/(twin-t0)), PFCp_ref< (alfa/(twin-t0)) | PFCp_tar< (alfa/(twin-t0))]

end



if plot_ctx_figures
  %% plot 3: rCLTC vs tCLT comparison with Passive
  
  
  figure;
  set(gcf, 'Position', [400 52 509 739]);
  
  
  % A1
  [mA1rr seA1rr] = jackmeanerr(A1_rCLT(A1tpr, :), 20);
  [mA1tt seA1tt] = jackmeanerr(A1_tCLT(A1tpr, :), 20);
  
  subplot(4,2,1);
  shadedErrorBar(A1_t, mA1rr, seA1rr, 'b', 0);
  hold on;
  shadedErrorBar(A1_t, mA1tt, seA1tt, 'r',   0);
  hold off;
  ylim(yrange);
  xlim([-0.4 1.8]);
  title('A1: tCLT vs rCLTC');
  
  subplot(4,2,2);
  shadedErrorBar(A1_t, mA1rp, seA1rp, 'b', 0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold on;
  shadedErrorBar(A1_t, mA1tp, seA1tp, 'r',   0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold off;
  xlim([-0.4 1.8]);
  title('A1: Passive ref vs tar');
  text(1.2, 0.5, ['N=' num2str(length(A1tpr))]);
  
  
  % dPEG
  if svd_subset == 0
    [mdPEGrr sedPEGrr] = jackmeanerr(dPEG_rCLT(dPEGtpr, :), 20);
    [mdPEGtt sedPEGtt] = jackmeanerr(dPEG_tCLT(dPEGtpr, :), 20);
    
    subplot(4,2,3);
    shadedErrorBar(dPEG_t, mdPEGrr, sedPEGrr, 'b', 0);
    hold on;
    shadedErrorBar(dPEG_t, mdPEGtt, sedPEGtt, 'r',   0);
    hold off;
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    xlim([-0.4 1.8]);
    title('dPEG: tCLT vs rCLTC');
    
    subplot(4,2,4);
    shadedErrorBar(dPEG_t, mdPEGrp, sedPEGrp, 'b', 0);
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    hold on;
    shadedErrorBar(dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    hold off;
    xlim([-0.4 1.8]);
    title('dPEG: Passive ref vs tar');
    text(1.2, 0.5, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % svd_dPEG
  if svd_subset == 1
    [mdPEGrr sedPEGrr] = jackmeanerr(svd_dPEG_rCLT(dPEGtpr, :), 20);
    [mdPEGtt sedPEGtt] = jackmeanerr(svd_dPEG_tCLT(dPEGtpr, :), 20);
    
    subplot(4,2,3);
    shadedErrorBar(svd_dPEG_t, mdPEGrr, sedPEGrr, 'b', 0);
    %ylim([-0.2 0.6]);
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtt, sedPEGtt, 'r',   0);
    hold off;
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    xlim([-0.4 1.8]);
    title('dPEG: tCLT vs rCLTC');
    
    subplot(4,2,4);
    shadedErrorBar(svd_dPEG_t, mdPEGrp, sedPEGrp, 'b', 0);
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    hold off;
    xlim([-0.4 1.8]);
    title('dPEG: Passive ref vs tar');
    text(1.2, 0.5, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % VPr
  [mVPrrr seVPrrr] = jackmeanerr(VPr_rCLT(VPrtpr, :), 20);
  [mVPrtt seVPrtt] = jackmeanerr(VPr_tCLT(VPrtpr, :), 20);
  
  subplot(4,2,5);
  shadedErrorBar(VPr_t, mVPrrr, seVPrrr, 'b', 0);
  hold on;
  shadedErrorBar(VPr_t, mVPrtt, seVPrtt, 'r',   0);
  hold off;
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  xlim([-0.4 1.8]);
  title('VPr: tCLT vs rCLTC');
  
  subplot(4,2,6);
  shadedErrorBar(VPr_t, mVPrrp, seVPrrp, 'b', 0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold on;
  shadedErrorBar(VPr_t, mVPrtp, seVPrtp, 'r',   0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold off;
  xlim([-0.4 1.8]);
  title('VPr: Passive ref vs tar');
  text(0.1, 0.5, ['N=' num2str(length(VPrtpr))]);
  
  
  % PFC
  [mPFCrr sePFCrr] = jackmeanerr(PFC_rCLT(PFCtpr, :), 20);
  [mPFCtt sePFCtt] = jackmeanerr(PFC_tCLT(PFCtpr, :), 20);
  
  subplot(4,2,7);
  shadedErrorBar(PFC_t, mPFCrr, sePFCrr, 'b', 0);
  hold on;
  shadedErrorBar(PFC_t, mPFCtt, sePFCtt, 'r',   0);
  hold off;
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  xlim([-0.4 1.8]);
  title('PFC: tCLT vs rCLTC');
  
  subplot(4,2,8);
  shadedErrorBar(PFC_t, mPFCrp, sePFCrp, 'b', 0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold on;
  shadedErrorBar(PFC_t, mPFCtp, sePFCtp, 'r',   0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  xlim([-0.4 1.8]);
  hold off;
  title('PFC: Passive ref vs tar');
  text(0.1, 0.5, ['N=' num2str(length(PFCtpr))]);
  
  suptitle('rCLTC vs. tone, out and in task context');
  
  set(gcf, 'PaperPositionMode', 'auto');
  print('03', '-dpdf');
  %% plot 4: rCLTC vs pass ref / tCLT vs pass tar
  
  
  figure;
  set(gcf, 'Position', [500 52 509 739]);
  
  
  % A1
  subplot(4,2,1);
  shadedErrorBar(A1_t, mA1rr, seA1rr, '--b', 0);
  hold on;
  shadedErrorBar(A1_t, mA1rp, seA1rp, 'b',   0);
  hold off;
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  xlim([-0.4 1.8]);
  title('A1: rCLTC vs pass rCLTC');
  
  subplot(4,2,2);
  shadedErrorBar(A1_t, mA1tt, seA1tt, '--r', 0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold on;
  shadedErrorBar(A1_t, mA1tp, seA1tp, 'r',   0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold off;
  xlim([-0.4 1.8]);
  title('A1: tCLT vs pass target');
  text(0.1, 0.5, ['N=' num2str(length(A1tpr))]);
  
  
  % dPEG
  if svd_subset == 0
    subplot(4,2,3);
    shadedErrorBar(dPEG_t, mdPEGrr, sedPEGrr, '--b', 0);
    hold on;
    shadedErrorBar(dPEG_t, mdPEGrp, sedPEGrp, 'b',   0);
    hold off;
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    xlim([-0.4 1.8]);
    title('dPEG: rCLTC vs pass rCLTC');
    
    subplot(4,2,4);
    shadedErrorBar(dPEG_t, mdPEGtt, sedPEGtt, '--r', 0);
    hold on;
    shadedErrorBar(dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    hold off;
    xlim([-0.4 1.8]);
    title('dPEG: tCLT vs pass target');
    text(0.1, 0.5, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % svd_dPEG
  if svd_subset == 1
    subplot(4,2,3);
    shadedErrorBar(svd_dPEG_t, mdPEGrr, sedPEGrr, '--b', 0);
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, 'b',   0);
    hold off;
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    xlim([-0.4 1.8]);
    title('dPEG: rCLTC vs pass rCLTC');
    
    subplot(4,2,4);
    shadedErrorBar(svd_dPEG_t, mdPEGtt, sedPEGtt, '--r', 0);
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    if only_pos_tar == -1
      ylim([-0.5 0.5]);
    else
      ylim([-0.2 0.6]);
    end
    hold off;
    xlim([-0.4 1.8]);
    title('dPEG: tCLT vs pass target');
    text(0.1, 0.5, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % VPr
  subplot(4,2,5);
  shadedErrorBar(VPr_t, mVPrrr, seVPrrr, '--b', 0);
  hold on;
  shadedErrorBar(VPr_t, mVPrrp, seVPrrp, 'b',   0);
  hold off;
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  xlim([-0.4 1.8]);
  title('VPr: rCLTC vs pass rCLTC');
  
  subplot(4,2,6);
  shadedErrorBar(VPr_t, mVPrtt, seVPrtt, '--r', 0);
  hold on;
  shadedErrorBar(VPr_t, mVPrtp, seVPrtp, 'r',   0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  hold off;
  xlim([-0.4 1.8]);
  title('VPr: tCLT vs pass target');
  text(0.1, 0.5, ['N=' num2str(length(VPrtpr))]);
  
  
  % PFC
  subplot(4,2,7);
  shadedErrorBar(PFC_t, mPFCrr, sePFCrr, '--b', 0);
  hold on;
  shadedErrorBar(PFC_t, mPFCrp, sePFCrp, 'b',   0);
  hold off;
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  xlim([-0.4 1.8]);
  title('PFC: rCLTC vs pass rCLTC');
  
  subplot(4,2,8);
  shadedErrorBar(PFC_t, mPFCtt, sePFCtt, '--r', 0);
  hold on;
  shadedErrorBar(PFC_t, mPFCtp, sePFCtp, 'r',   0);
  if only_pos_tar == -1
    ylim([-0.5 0.5]);
  else
    ylim([-0.2 0.6]);
  end
  xlim([-0.4 1.8]);
  hold off;
  title('PFC: tCLT vs pass target');
  text(0.1, 0.5, ['N=' num2str(length(PFCtpr))]);
  
  suptitle('Out of task vs. passive task');
  
  
  set(gcf, 'PaperPositionMode', 'auto');
  print('04', '-dpdf');
end

%% Plot 5: Naive vs Trained figure
if add_naive == 1
  
  if only_pos_tar == -1
    yrange([-0.5 0.5]);
  else
    yrange =[-0.2 0.6];
  end
  
  if popnorm
    yrange = [-1 1.3];
  elseif ~popnorm && ~normalize_each
    yrange = [-20 20];
  end
  
  figure;
  set(gcf, 'Position', [21 48 693 741]);
  set(gcf, 'Color', 'w');
  % A1

  subplot(3,3,1);
  h5 = shadedErrorBar(nA1_t, mnA1rp, senA1rp, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(nA1_t, mnA1tp, senA1tp, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  ylim(yrange);
  xlim([0.9 2.2]);
  hold off;
  title('\bfNaive Passive CRD');
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(nA1_cells_used)))]);
  set(gca, 'box', 'off');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfA1', 'FontSize', 12);
  set(gca, 'XTickLabel', []);
  set(gca, 'Layer', 'top');

  subplot(3,3,2);
  h5 = shadedErrorBar(A1_t, mA1rp, seA1rp, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(A1_t, mA1tp, seA1tp, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  title('\bfTrained Passive CRD');
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(A1tpr))]);
  
  subplot(3,3,3);
  h5 = shadedErrorBar(A1_t, mA1ra, seA1ra, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(A1_t, mA1ta, seA1ta, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([1 2.2]);
  ylim(yrange);
  title('\bfBehaving CRD');
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  
  % dPEG
  if svd_subset == 1
    t = svd_dPEG_t;
  else
    t = dPEG_t;
  end
  

  subplot(3,3,4);
  h5 = shadedErrorBar(ndPEG_t, mndPEGrp, sendPEGrp, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(ndPEG_t, mndPEGtp, sendPEGtp, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold off;
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  ylim(yrange);
  xlim([.9 2.2]);
  hold off;
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(ndPEG_cells_used)))]);
  set(gca, 'box', 'off');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfdPEG', 'FontSize', 12);
  set(gca, 'XTickLabel', []);
  set(gca, 'Layer', 'top');

  subplot(3,3,5);
  h5 = shadedErrorBar(t, mdPEGrp, sedPEGrp, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(t, mdPEGtp, sedPEGtp, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(dPEGtpr))]);
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,6);
  h5 = shadedErrorBar(t, mdPEGra, sedPEGra, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(t, mdPEGta, sedPEGta, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([1 2.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  
  % VPr
  subplot(3,3,7);
  h5 = shadedErrorBar(nVPr_t, mnVPrrp, senVPrrp, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(nVPr_t, mnVPrtp, senVPrtp, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  xlim([.9 2.2]);
  ylim(yrange);
  hold off;
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(nVPr_cells_used)))]);
  set(gca, 'box', 'off');
  set(gca, 'Layer', 'top');
  xlabel('Time (s)');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfVPr', 'FontSize', 12);

  subplot(3,3,8);
  h5 = shadedErrorBar(VPr_t, mVPrrp, seVPrrp, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(VPr_t, mVPrtp, seVPrtp, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(VPrtpr))]);
  set(gca, 'box', 'off');
  xlabel('Time (s)');
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');

  subplot(3,3,9);
  h5 = shadedErrorBar(VPr_t, mVPrra, seVPrra, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  hold on;
  h5 = shadedErrorBar(VPr_t, mVPrta, seVPrta, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  set(gca, 'box', 'off');
  xlabel('Time (s)');
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  
  set(gcf, 'PaperPositionMode', 'auto');
  print('05', '-dpdf');
end

%% Plot 6: Tar-Ref differences

% A1

% dif between target and reference in passive and active
tempdiff1 = A1_Tpassive(A1tpr,:) - A1_Rpassive(A1tpr,:);
[A1_Dtr_passive  A1_Dtr_passive_sem] = jackmeanerr(tempdiff1, 20); % delta tar-ref
A1_Dtrp = nanmean(tempdiff1(:,scatwin),2);

tempdiff2 = A1_Tactive(A1tpr,:) - A1_Ractive(A1tpr,:);
[A1_Dtr_active A1_Dtr_active_sem]    = jackmeanerr(tempdiff2, 20);
A1_Dtra = nanmean(tempdiff2(:,scatwin),2);

A1_contrastenh = tempdiff2 - tempdiff1; % difference of the difference; enhancement of the tar-ref contrast in behavior
[mA1_contrastenh seA1_contrastenh] = jackmeanerr(A1_contrastenh, 20);

if add_naive == 1
  ntempdiff = nA1_Tpassive - nA1_Rpassive; % naive
  [nA1_Dtr_passive  nA1_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
  nA1_Dtrp = nanmean(ntempdiff(:,nscatwin),2);
end

% diff between active and passive in reference and target
tempdiff1 = A1_Ractive(A1tpr,:) - A1_Rpassive(A1tpr,:);
[A1_Dap_ref A1_Dap_ref_sem]    = jackmeanerr(tempdiff1, 20); 
A1_Dapr = nanmean(tempdiff1(:,scatwin),2);

tempdiff2 = A1_Tactive(A1tpr,:) - A1_Tpassive(A1tpr,:);
[A1_Dap_tar A1_Dap_tar_sem]    = jackmeanerr(tempdiff2, 20); 
A1_Dapt = nanmean(tempdiff2(:,scatwin),2);
A1_at = nanmean(A1_Tactive(A1tpr, scatwin),2);

A1_tarenhancement = tempdiff2 - tempdiff1; % difference of the difference; target change in active minus reference change
[mA1_tarenhancement seA1_tarenhancement] = jackmeanerr(A1_tarenhancement, 20);

% dPEG
if svd_subset == 0
  % dif between target and reference in passive and active
  tempdiff1 = dPEG_Tpassive(dPEGtpr,:) - dPEG_Rpassive(dPEGtpr,:);
  [dPEG_Dtr_passive  dPEG_Dtr_passive_sem] = jackmeanerr(tempdiff1, 20); % delta tar-ref
  dPEG_Dtrp = nanmean(tempdiff1(:,scatwin),2);
  
  tempdiff2 = dPEG_Tactive(dPEGtpr,:) - dPEG_Ractive(dPEGtpr,:);
  [dPEG_Dtr_active dPEG_Dtr_active_sem]    = jackmeanerr(tempdiff2, 20);
  dPEG_Dtra = nanmean(tempdiff2(:,scatwin),2);
  
  dPEG_contrastenh = tempdiff2 - tempdiff1; % difference of the difference; enhancement of the tar-ref contrast in behavior
  [mdPEG_contrastenh sedPEG_contrastenh] = jackmeanerr(dPEG_contrastenh, 20);
  
  if add_naive == 1
    ntempdiff = ndPEG_Tpassive - ndPEG_Rpassive; % naive
    [ndPEG_Dtr_passive  ndPEG_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
    ndPEG_Dtrp = nanmean(ntempdiff(:,nscatwin),2);
  end
  
  % diff between active and passive in reference and target
  tempdiff1 = dPEG_Ractive(dPEGtpr,:) - dPEG_Rpassive(dPEGtpr,:);
  [dPEG_Dap_ref dPEG_Dap_ref_sem]    = jackmeanerr(tempdiff1, 20);
  dPEG_Dapr = nanmean(tempdiff1(:,scatwin),2);
  
  tempdiff2 = dPEG_Tactive(dPEGtpr,:) - dPEG_Tpassive(dPEGtpr,:);
  [dPEG_Dap_tar dPEG_Dap_tar_sem]    = jackmeanerr(tempdiff2, 20);
  dPEG_Dapt = nanmean(tempdiff2(:,scatwin),2);
  dPEG_at = nanmean(dPEG_Tactive(dPEGtpr, scatwin),2);
  
  dPEG_tarenhancement = tempdiff2 - tempdiff1; % difference of the difference; target change in active minus reference change
  [mdPEG_tarenhancement sedPEG_tarenhancement] = jackmeanerr(dPEG_tarenhancement, 20);
  
elseif svd_subset == 1
  % svd_dPEG
  % dif between target and reference in passive and active
  tempdiff1 = svd_dPEG_Tpassive(dPEGtpr,:) - svd_dPEG_Rpassive(dPEGtpr,:);
  [dPEG_Dtr_passive  dPEG_Dtr_passive_sem] = jackmeanerr(tempdiff1, 20); % delta tar-ref
  dPEG_Dtrp = nanmean(tempdiff1(:,scatwin),2);
  
  tempdiff2 = svd_dPEG_Tactive(dPEGtpr,:) - svd_dPEG_Ractive(dPEGtpr,:);
  [dPEG_Dtr_active dPEG_Dtr_active_sem]    = jackmeanerr(tempdiff2, 20);
  dPEG_Dtra = nanmean(tempdiff2(:,scatwin),2);
  
  dPEG_contrastenh = tempdiff2 - tempdiff1; % difference of the difference; enhancement of the tar-ref contrast in behavior
  [mdPEG_contrastenh sedPEG_contrastenh] = jackmeanerr(dPEG_contrastenh, 20);
  
  if add_naive == 1
    ntempdiff = ndPEG_Tpassive - ndPEG_Rpassive; % naive
    [ndPEG_Dtr_passive  ndPEG_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
    ndPEG_Dtrp = nanmean(ntempdiff(:,nscatwin),2);
  end
  
  % diff between active and passive in reference and target
  tempdiff1 = svd_dPEG_Ractive(dPEGtpr,:) - svd_dPEG_Rpassive(dPEGtpr,:);
  [dPEG_Dap_ref dPEG_Dap_ref_sem]    = jackmeanerr(tempdiff1, 20);
  dPEG_Dapr = nanmean(tempdiff1(:,scatwin),2);
  
  tempdiff2 = svd_dPEG_Tactive(dPEGtpr,:) - svd_dPEG_Tpassive(dPEGtpr,:);
  [dPEG_Dap_tar dPEG_Dap_tar_sem]    = jackmeanerr(tempdiff2, 20);
  dPEG_Dapt = nanmean(tempdiff2(:,scatwin),2);
  dPEG_at = nanmean(dPEG_Tactive(dPEGtpr, scatwin),2);
  
  dPEG_tarenhancement = tempdiff2 - tempdiff1; % difference of the difference; target change in active minus reference change
  [mdPEG_tarenhancement sedPEG_tarenhancement] = jackmeanerr(dPEG_tarenhancement, 20);
  
end
  

% VPr
% dif between target and reference in passive and active
tempdiff1 = VPr_Tpassive(VPrtpr,:) - VPr_Rpassive(VPrtpr,:);
[VPr_Dtr_passive  VPr_Dtr_passive_sem] = jackmeanerr(tempdiff1, 20); % delta tar-ref
VPr_Dtrp = nanmean(tempdiff1(:,scatwin),2);

tempdiff2 = VPr_Tactive(VPrtpr,:) - VPr_Ractive(VPrtpr,:);
[VPr_Dtr_active VPr_Dtr_active_sem]    = jackmeanerr(tempdiff2, 20);
VPr_Dtra = nanmean(tempdiff2(:,scatwin),2);


VPr_contrastenh = tempdiff2 - tempdiff1; % difference of the difference; enhancement of the tar-ref contrast in behavior
[mVPr_contrastenh seVPr_contrastenh] = jackmeanerr(VPr_contrastenh, 20);

if add_naive == 1
  ntempdiff = nVPr_Tpassive - nVPr_Rpassive; % naive
  [nVPr_Dtr_passive  nVPr_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
  nVPr_Dtrp = nanmean(ntempdiff(:,nscatwin),2);
end

% diff between active and passive in reference and target
tempdiff1 = VPr_Ractive(VPrtpr,:) - VPr_Rpassive(VPrtpr,:);
[VPr_Dap_ref VPr_Dap_ref_sem]    = jackmeanerr(tempdiff1, 20); 
VPr_Dapr = nanmean(tempdiff1(:,scatwin),2);

tempdiff2 = VPr_Tactive(VPrtpr,:) - VPr_Tpassive(VPrtpr,:);
[VPr_Dap_tar VPr_Dap_tar_sem]    = jackmeanerr(tempdiff2, 20); 
VPr_Dapt = nanmean(tempdiff2(:,scatwin),2);
VPr_at = nanmean(VPr_Tactive(VPrtpr, scatwin),2);

VPr_tarenhancement = tempdiff2 - tempdiff1; % difference of the difference; target change in active minus reference change
[mVPr_tarenhancement seVPr_tarenhancement] = jackmeanerr(VPr_tarenhancement, 20);


% PFC

% dif between target and reference in passive and active
tempdiff1 = PFC_Tpassive(PFCtpr,:) - PFC_Rpassive(PFCtpr,:);
[PFC_Dtr_passive  PFC_Dtr_passive_sem] = jackmeanerr(tempdiff1, 20); % delta tar-ref
PFC_Dtrp = nanmean(tempdiff1(:,scatwin),2);

tempdiff2 = PFC_Tactive(PFCtpr,:) - PFC_Ractive(PFCtpr,:);
[PFC_Dtr_active PFC_Dtr_active_sem]    = jackmeanerr(tempdiff2, 20);
PFC_Dtra = nanmean(tempdiff2(:,scatwin),2);

PFC_contrastenh = tempdiff2 - tempdiff1; % difference of the difference; enhancement of the tar-ref contrast in behavior
[mPFC_contrastenh sePFC_contrastenh] = jackmeanerr(PFC_contrastenh, 20);

% diff between active and passive in reference and target
tempdiff1 = PFC_Ractive(PFCtpr,:) - PFC_Rpassive(PFCtpr,:);
[PFC_Dap_ref PFC_Dap_ref_sem]    = jackmeanerr(tempdiff1, 20); 
PFC_Dapr = nanmean(tempdiff1(:,scatwin),2);

tempdiff2 = PFC_Tactive(PFCtpr,:) - PFC_Tpassive(PFCtpr,:);
[PFC_Dap_tar PFC_Dap_tar_sem]    = jackmeanerr(tempdiff2, 20); 
PFC_Dapt = nanmean(tempdiff2(:,scatwin),2);
PFC_at = nanmean(PFC_Tactive(PFCtpr, scatwin),2);

PFC_tarenhancement = tempdiff2 - tempdiff1; % difference of the difference; target change in active minus reference change
[mPFC_tarenhancement sePFC_tarenhancement] = jackmeanerr(PFC_tarenhancement, 20);

% definitions done, plotting

figure;
set(gcf, 'Position', [21 48 693 741]);

% A1
subplot(4,3,1);
shadedErrorBar(A1_t, A1_Dap_ref, A1_Dap_ref_sem, 'b', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1 Change in reference response');

subplot(4,3,2);
shadedErrorBar(A1_t, A1_Dap_tar, A1_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Change in target response');

subplot(4,3,3);
shadedErrorBar(A1_t, mA1_tarenhancement, seA1_tarenhancement, 'k', 0);
xlim([-0.4 3.2]);
title('Tar change - Ref change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(A1tpr))]);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);


% dPEG
subplot(4,3,4);
shadedErrorBar(dPEG_t, dPEG_Dap_ref, dPEG_Dap_ref_sem, 'b', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG Change in reference response');

subplot(4,3,5);
shadedErrorBar(dPEG_t, dPEG_Dap_tar, dPEG_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Change in target response');

subplot(4,3,6);
shadedErrorBar(dPEG_t, mdPEG_tarenhancement, sedPEG_tarenhancement, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Tar change - Ref change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);


% VPr
subplot(4,3,7);
shadedErrorBar(VPr_t, VPr_Dap_ref, VPr_Dap_ref_sem, 'b', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr Change in reference response');

subplot(4,3,8);
shadedErrorBar(VPr_t, VPr_Dap_tar, VPr_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Change in target response');

subplot(4,3,9);
shadedErrorBar(VPr_t, mVPr_tarenhancement, seVPr_tarenhancement, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Tar change - Ref change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(VPrtpr))]);


% PFC
subplot(4,3,10);
shadedErrorBar(PFC_t, PFC_Dap_ref, PFC_Dap_ref_sem, 'b', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off
xlim([-0.4 3.2]);
ylim(yrange);
title('PFC Change in reference response');

subplot(4,3,11);
shadedErrorBar(PFC_t, PFC_Dap_tar, PFC_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Change in target response');

subplot(4,3,12);
shadedErrorBar(PFC_t, mPFC_tarenhancement, sePFC_tarenhancement, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Tar change - Ref change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(PFCtpr))]);

suptitle('Active - Passive');

set(gcf, 'PaperPositionMode', 'auto');
print('06', '-dpdf');

%% figure 7: target minus reference area comparison

figure;
set(gcf, 'Position', [21 48 693 741]);

% A1
subplot(4,3,1);
shadedErrorBar(A1_t, A1_Dtr_passive, A1_Dtr_passive_sem, '--m', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1 Passive tar-ref contrast');

subplot(4,3,2);
shadedErrorBar(A1_t, A1_Dtr_active, A1_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Active tar-ref contrast');

subplot(4,3,3);
shadedErrorBar(A1_t, mA1_contrastenh, seA1_contrastenh, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Tar - Ref contrast change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(A1tpr))]);


% dPEG
subplot(4,3,4);
shadedErrorBar(dPEG_t, dPEG_Dtr_passive, dPEG_Dtr_passive_sem, '--m', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG Passive tar-ref contrast');

subplot(4,3,5);
shadedErrorBar(dPEG_t, dPEG_Dtr_active, dPEG_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Active tar-ref contrast');

subplot(4,3,6);
shadedErrorBar(dPEG_t, mdPEG_contrastenh, sedPEG_contrastenh, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Tar - Ref contrast change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);


% VPr
subplot(4,3,7);
shadedErrorBar(VPr_t, VPr_Dtr_passive, VPr_Dtr_passive_sem, '--m', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr Passive tar-ref contrast');

subplot(4,3,8);
shadedErrorBar(VPr_t, VPr_Dtr_active, VPr_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Active tar-ref contrast');

subplot(4,3,9);
shadedErrorBar(VPr_t, mVPr_contrastenh, seVPr_contrastenh, 'k', 0);
xlim([-0.4 3.2]);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8], yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Tar - Ref contrast change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(VPrtpr))]);


% PFC
subplot(4,3,10);
shadedErrorBar(PFC_t, PFC_Dtr_passive, PFC_Dtr_passive_sem, '--m', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('PFC Passive tar-ref contrast');

subplot(4,3,11);
shadedErrorBar(PFC_t, PFC_Dtr_active, PFC_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Active tar-ref contrast');

subplot(4,3,12);
shadedErrorBar(PFC_t, mPFC_contrastenh, sePFC_contrastenh, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('Tar - Ref contrast change');
text(0.1, yrange(2)*.8, ['N=' num2str(length(PFCtpr))]);

suptitle('Tar-Ref contrast in passive and active');

set(gcf, 'PaperPositionMode', 'auto');
print('07', '-dpdf');

%% figure 8: scatter plots ala fig 5 in neuron 2014
if only_pos_tar == -1
  yrange([-1 1]);
else
  yrange =[-1 1]; % i know this is dumb...
end

if popnorm
  yrange = [-5 5];
elseif popnorm == 0 && normalize_each == 0
  yrange = [-15 15];
  
end

figure;
subplot(2,2,1);
scatter(A1_Dtrp, A1_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange);
ylim(yrange);
xlabel('Passive target-reference');
ylabel('Active target-reference');
title('A1');
hold on;
line(yrange, yrange, 'Color', 'k');
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k');
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k');

subplot(2,2,2);
scatter(dPEG_Dtrp, dPEG_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange);
ylim(yrange);
xlabel('Passive target-reference');
ylabel('Active target-reference');
title('dPEG');
hold on
line(yrange, yrange, 'Color', 'k');
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k');
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k');

subplot(2,2,3);
scatter(VPr_Dtrp, VPr_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange);
ylim(yrange);
xlabel('Passive target-reference');
ylabel('Active target-reference');
title('VPr');
hold on;
line(yrange, yrange, 'Color', 'k');
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k');
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k');

subplot(2,2,4);
scatter(PFC_Dtrp, PFC_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange);
ylim(yrange);
xlabel('Passive target-reference');
ylabel('Active target-reference');
title('PFC');
hold on;
line(yrange, yrange, 'Color', 'k');
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k');
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k');

set(gcf, 'PaperPositionMode', 'auto');
print('09', '-dpdf');

%% figure 9 : target enhancement index -- histogram
edges = [-1:0.1:1];

if normalize_each == 0
  edges = [-10: 2: 10];
end
if popnorm == 1
  edges = [-3:0.25:3];
  xrange = [-3 3];
elseif popnorm == 0 && normalize_each == 0
  xrange = [-20 20];
  edges = [-20: 2: 20];
elseif popnorm == 0 && normalize_each == 1
  xrange = [-1 1];
  edges = [-1: .2: 1];
end

figure;
subplot(2,2,1);
A1_h = histc(A1_Dtra - A1_Dtrp, edges);
hold on;
bar(edges, A1_h, 'style', 'histc');
line([0 0], [0 max(A1_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
hold off;
xlim(xrange);
ylim([0 max(A1_h)*1.2]);
ylabel('Number of cells');
title('A1');

subplot(2,2,2);
dPEG_h = histc(dPEG_Dtra - dPEG_Dtrp, edges);
hold on;
bar(edges, dPEG_h, 'style', 'histc');
line([0 0], [0 max(dPEG_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
hold off;
xlim(xrange);
ylim([0 max(dPEG_h)*1.2]);
ylabel('Number of cells');
title('dPEG');

subplot(2,2,3);
VPr_h = histc(VPr_Dtra - VPr_Dtrp, edges);
hold on;
bar(edges, VPr_h, 'style', 'histc');
line([0 0], [0 max(VPr_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
hold off;
xlim(xrange);
ylim([0 max(VPr_h)*1.2]);
ylabel('Number of cells');
title('VPr');

subplot(2,2,4);
PFC_h = histc(PFC_Dtra - PFC_Dtrp, edges);
hold on;
bar(edges, PFC_h, 'style', 'histc');
line([0 0], [0 max(PFC_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
hold off;
xlim(xrange);
ylim([0 max(PFC_h)*1.2]);
ylabel('Number of cells');
title('dlFC');

set(gcf, 'PaperPositionMode', 'auto');
print('10', '-dpdf');

%% figure 10 : naive, passive and active tar-ref contrast histograms

lefthlf = [1: floor(length(edges)/2)];
righthlf = [floor(length(edges)/2)+1 : length(edges)];
leftclr = [190/255, 65/255, 190/255];
rightclr = [128/255,0/255,128/255];

figure;
set(gcf, 'Position', [2 77 559 747]);
set(gcf, 'Color', 'w');

subplot(4,3,1);
nA1_h = histc(nA1_Dtrp, edges);
hold on;
ngb = bar(edges(lefthlf), nA1_h(lefthlf), 'style', 'histc');
set(ngb, 'FaceColor', leftclr);
set(ngb, 'EdgeColor', leftclr);
nbb = bar(edges(righthlf), nA1_h(righthlf), 'style', 'histc');
set(nbb, 'FaceColor', rightclr);
set(nbb, 'EdgeColor', rightclr);
line([0 0], [0 max(nA1_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(nA1_Dtrp), max(nA1_h)*1.1, num2str(nanmean(nA1_Dtrp),1), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(nA1_h)*1.2]);
ylabel('Number of cells');
text(-6, max(nA1_h)*.6, '\bfA1', 'FontSize', 11);
title('\bfNaive');
set(gca, 'Layer', 'top');

subplot(4,3,2);
A1_h = histc(A1_Dtrp, edges);
hold on;
pgb = bar(edges(lefthlf), A1_h(lefthlf), 'style', 'histc');
set(pgb, 'FaceColor', leftclr);
set(pgb, 'EdgeColor', leftclr);
pbb = bar(edges(righthlf), A1_h(righthlf), 'style', 'histc');
set(pbb, 'FaceColor', rightclr);
set(pbb, 'EdgeColor', rightclr);
line([0 0], [0 max(A1_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(A1_Dtrp), max(A1_h)*1.1, num2str(nanmean(A1_Dtrp),1), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(A1_h)*1.2]);
title('\bfTrained Passive');
set(gca, 'Layer', 'top');

subplot(4,3,3);
A1_h = histc(A1_Dtra, edges);
hold on;
agb = bar(edges(lefthlf), A1_h(lefthlf), 'style', 'histc');
set(agb, 'FaceColor', leftclr);
set(agb, 'EdgeColor', leftclr);
abb = bar(edges(righthlf), A1_h(righthlf), 'style', 'histc');
set(abb, 'FaceColor', rightclr);
set(abb, 'EdgeColor', rightclr);
line([0 0], [0 max(A1_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(A1_Dtra), max(A1_h)*1.1, num2str(nanmean(A1_Dtra),1), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(A1_h)*1.2]);
title('\bfBehaving');
set(gca, 'Layer', 'top');

subplot(4,3,4);
ndPEG_h = histc(ndPEG_Dtrp, edges);
hold on;
ngb = bar(edges(lefthlf), ndPEG_h(lefthlf), 'style', 'histc');
set(ngb, 'FaceColor', leftclr);
set(ngb, 'EdgeColor', leftclr);
nbb = bar(edges(righthlf), ndPEG_h(righthlf), 'style', 'histc');
set(nbb, 'FaceColor', rightclr);
set(nbb, 'EdgeColor', rightclr);
line([0 0], [0 max(ndPEG_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(ndPEG_Dtrp), max(ndPEG_h)*1.1, num2str(nanmean(ndPEG_Dtrp),1), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(ndPEG_h)*1.2]);
ylabel('Number of cells');
text(-7, max(ndPEG_h)*.6, '\bfdPEG', 'FontSize', 11);
set(gca, 'Layer', 'top');

subplot(4,3,5);
dPEG_h = histc(dPEG_Dtrp, edges);
hold on;
pgb = bar(edges(lefthlf), dPEG_h(lefthlf), 'style', 'histc');
set(pgb, 'FaceColor', leftclr);
set(pgb, 'EdgeColor', leftclr);
pbb = bar(edges(righthlf), dPEG_h(righthlf), 'style', 'histc');
set(pbb, 'FaceColor', rightclr);
set(pbb, 'EdgeColor', rightclr);
line([0 0], [0 max(dPEG_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(dPEG_Dtrp), max(dPEG_h)*1.1, num2str(nanmean(dPEG_Dtrp),2), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(dPEG_h)*1.2]);
set(gca, 'Layer', 'top');

subplot(4,3,6);
dPEG_h = histc(dPEG_Dtra, edges);
hold on;
agb = bar(edges(lefthlf), dPEG_h(lefthlf), 'style', 'histc');
set(agb, 'FaceColor', leftclr);
set(agb, 'EdgeColor', leftclr);
abb = bar(edges(righthlf), dPEG_h(righthlf), 'style', 'histc');
set(abb, 'FaceColor', rightclr);
set(abb, 'EdgeColor', rightclr);
line([0 0], [0 max(dPEG_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(dPEG_Dtra), max(dPEG_h)*1.1, num2str(nanmean(dPEG_Dtra),2), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(dPEG_h)*1.2]);
set(gca, 'Layer', 'top');

subplot(4,3,7);
nVPr_h = histc(nVPr_Dtrp, edges);
hold on;
ngb = bar(edges(lefthlf), nVPr_h(lefthlf), 'style', 'histc');
set(ngb, 'FaceColor', leftclr);
set(ngb, 'EdgeColor', leftclr);
nbb = bar(edges(righthlf), nVPr_h(righthlf), 'style', 'histc');
set(nbb, 'FaceColor', rightclr);
set(nbb, 'EdgeColor', rightclr);
line([0 0], [0 max(nVPr_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(nVPr_Dtrp), max(nVPr_h)*1.1, num2str(nanmean(nVPr_Dtrp),2), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(nVPr_h)*1.2]);
ylabel('Number of cells');
xlabel(sprintf('Fraction Target/Reference\ncontrast'));
text(-6.5, max(nVPr_h)*.6, '\bfVPr', 'FontSize', 11);
set(gca, 'Layer', 'top');

subplot(4,3,8);
VPr_h = histc(VPr_Dtrp, edges);
hold on;
pgb = bar(edges(lefthlf), VPr_h(lefthlf), 'style', 'histc');
set(pgb, 'FaceColor', leftclr);
set(pgb, 'EdgeColor', leftclr);
pbb = bar(edges(righthlf), VPr_h(righthlf), 'style', 'histc');
set(pbb, 'FaceColor', rightclr);
set(pbb, 'EdgeColor', rightclr);
line([0 0], [0 max(VPr_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(VPr_Dtrp), max(VPr_h)*1.1, num2str(nanmean(VPr_Dtrp),2), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(VPr_h)*1.2]);
set(gca, 'Layer', 'top');

subplot(4,3,9);
VPr_h = histc(VPr_Dtra, edges);
% pd = fitdist(VPr_Dtra, 'normal');
% fd = pdf(pd, edges);
% fd = fd* (max(VPr_h)*3);
hold on;
agb = bar(edges(lefthlf), VPr_h(lefthlf), 'style', 'histc');
set(agb, 'FaceColor', leftclr);
set(agb, 'EdgeColor', leftclr);
abb = bar(edges(righthlf), VPr_h(righthlf), 'style', 'histc');
set(abb, 'FaceColor', rightclr);
set(abb, 'EdgeColor', rightclr);
line([0 0], [0 max(VPr_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(VPr_Dtra), max(VPr_h)*1.1, num2str(nanmean(VPr_Dtra),2), 'Color', 'r');
% plot(edges,fd, 'r', 'LineWidth', 2);
hold off;
xlim(xrange);
ylim([0 max(VPr_h)*1.2]);
set(gca, 'Layer', 'top');

subplot(4,3,11);
PFC_h = histc(PFC_Dtrp, edges);
hold on;
pgb = bar(edges(lefthlf), PFC_h(lefthlf), 'style', 'histc');
set(pgb, 'FaceColor', leftclr);
set(pgb, 'EdgeColor', leftclr);
pbb = bar(edges(righthlf), PFC_h(righthlf), 'style', 'histc');
set(pbb, 'FaceColor', rightclr);
set(pbb, 'EdgeColor', rightclr);
line([0 0], [0 max(PFC_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(PFC_Dtrp), max(PFC_h)*1.1, num2str(nanmean(PFC_Dtrp),2), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(PFC_h)*1.2]);
ylabel('Number of cells');
xlabel(sprintf('Fraction Target/Reference\ncontrast'));
text(-15.0, max(PFC_h)*.6, '\bfdlFC', 'FontSize', 11);
set(gca, 'Layer', 'top');

subplot(4,3,12);
PFC_h = histc(PFC_Dtra, edges);
hold on;
agb = bar(edges(lefthlf), PFC_h(lefthlf), 'style', 'histc');
set(agb, 'FaceColor', leftclr);
set(agb, 'EdgeColor', leftclr);
abb = bar(edges(righthlf), PFC_h(righthlf), 'style', 'histc');
set(abb, 'FaceColor', rightclr);
set(abb, 'EdgeColor', rightclr);
line([0 0], [0 max(PFC_h)*1.2], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
text(nanmean(PFC_Dtra), max(PFC_h)*1.1, num2str(nanmean(PFC_Dtra),2), 'Color', 'r');
hold off;
xlim(xrange);
ylim([0 max(PFC_h)*1.2]);
xlabel(sprintf('Fraction Target/Reference\ncontrast'));
set(gca, 'Layer', 'top');


set(gcf, 'PaperPositionMode', 'auto');
print('11', '-dpdf');


%% Figure 11 definitions
%  Analysis of post-target response

% ptwin1 = [76:76+10];%89];%97]; %2.1 to 2.43 s from stim onset. 0.1 s after offset
% ptwin1 = [73:94]; %2.0 to 2.8 s from stim onset. 0 s after offset
ptwin1 = [74:97]; %2.05 to 2.8 s from stim onset. 0 s after offset
ptwin2 = ptwin1;

A1_ptp   = A1_Tpassive(A1tpr,ptwin1);
A1_pta   = A1_Tactive(A1tpr, ptwin1);
if svd_subset == 0
  dPEG_ptp = dPEG_Tpassive(dPEGtpr, ptwin1);
  dPEG_pta = dPEG_Tactive(dPEGtpr, ptwin1);
elseif svd_subset == 1
  dPEG_ptp = svd_dPEG_Tpassive(dPEGtpr, ptwin1);
  dPEG_pta = svd_dPEG_Tactive(dPEGtpr, ptwin1);
end
VPr_ptp  = VPr_Tpassive(VPrtpr, ptwin2);
VPr_pta  = VPr_Tactive(VPrtpr, ptwin2);
PFC_ptp  = PFC_Tpassive(PFCtpr, ptwin1);
PFC_pta  = PFC_Tactive(PFCtpr, ptwin1);

mA1_ptp = nanmean(A1_ptp, 2);
mA1_pta = nanmean(A1_pta, 2);
mdPEG_ptp = nanmean(dPEG_ptp, 2);
mdPEG_pta = nanmean(dPEG_pta, 2);
mVPr_ptp = nanmean(VPr_ptp, 2);
mVPr_pta = nanmean(VPr_pta, 2);
mPFC_ptp = nanmean(PFC_ptp, 2);
mPFC_pta = nanmean(PFC_pta, 2);

%post reference response
A1_prp   = A1_Rpassive(A1tpr,ptwin1);
A1_pra   = A1_Ractive(A1tpr, ptwin1);
if svd_subset == 0
  dPEG_prp = dPEG_Rpassive(dPEGtpr, ptwin1);
  dPEG_pra = dPEG_Ractive(dPEGtpr, ptwin1);
elseif svd_subset == 1
  dPEG_prp = svd_dPEG_Rpassive(dPEGtpr, ptwin1);
  dPEG_pra = svd_dPEG_Ractive(dPEGtpr, ptwin1);
end
VPr_prp  = VPr_Rpassive(VPrtpr, ptwin2);
VPr_pra  = VPr_Ractive(VPrtpr, ptwin2);
PFC_prp  = PFC_Rpassive(PFCtpr, ptwin1);
PFC_pra  = PFC_Ractive(PFCtpr, ptwin1);

mA1_prp = nanmean(A1_prp, 2);
mA1_pra = nanmean(A1_pra, 2);
mdPEG_prp = nanmean(dPEG_prp, 2);
mdPEG_pra = nanmean(dPEG_pra, 2);
mVPr_prp = nanmean(VPr_prp, 2);
mVPr_pra = nanmean(VPr_pra, 2);
mPFC_prp = nanmean(PFC_prp, 2);
mPFC_pra = nanmean(PFC_pra, 2);

[mA1_change seA1_change]     = jackmeanerr(mA1_pta-mA1_ptp);
[mdPEG_change sedPEG_change] = jackmeanerr(mdPEG_pta-mdPEG_ptp);
[mVPr_change seVPr_change]   = jackmeanerr(mVPr_pta-mVPr_ptp);
[mPFC_change sePFC_change]   = jackmeanerr(mPFC_pta-mPFC_ptp);

dA1_pt   = (mA1_pta-mA1_ptp);
ddPEG_pt = (mdPEG_pta-mdPEG_ptp);
dVPr_pt  = (mVPr_pta-mVPr_ptp);
dPFC_pt  = (mPFC_pta-mPFC_ptp);


%% figure 11
figure;
set(gcf, 'Position', [201 331 947 425]);
subplot(2,1,1);
hold on;
bar([mA1_change, mdPEG_change, mVPr_change, mPFC_change]);
errorbar([mA1_change, mdPEG_change, mVPr_change, mPFC_change], [seA1_change, sedPEG_change, seVPr_change, sePFC_change], '.r');
hold off;
set(gca, 'XTick', [1:4]);
set(gca, 'XTickLabel', {'A1', 'dPEG', 'VPr', 'dlFC'});
if normalize_each
  ylabel(sprintf('Normalized spike rate\nchange after target'));
else
  ylabel(sprintf('Spike rate\nchange after target'));
end



subplot(2,4,5);
hA1 = histc(mA1_pta-mA1_ptp, edges);
bar(edges, hA1, 'Style', 'hist');
hold on;
line([0 0], [0 max(hA1)*1.2], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
hold off;
xlim(xrange);
ylim([0 max(hA1)*1.2]);
title('A1');
ylabel('Number of cells');
xlabel(sprintf('Fraction change\nafter target'));

subplot(2,4,6);
hdPEG = histc(mdPEG_pta-mdPEG_ptp, edges);
bar(edges, hdPEG, 'Style', 'hist');
hold on;
line([0 0], [0 max(hdPEG)*1.2], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
hold off;
xlim(xrange);
ylim([0 max(hdPEG)*1.2]);
title('dPEG');
xlabel(sprintf('Fraction change\nafter target'));

subplot(2,4,7);
hVPr = histc(mVPr_pta-mVPr_ptp, edges);
bar(edges, hVPr, 'Style', 'hist');
hold on;
line([0 0], [0 max(hVPr)*1.2], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
hold off;
xlim(xrange);
ylim([0 max(hVPr)*1.2]);
title('VPr');
xlabel(sprintf('Fraction change\nafter target'));

subplot(2,4,8);
hPFC = histc(mPFC_pta-mPFC_ptp, edges);
bar(edges, hPFC, 'Style', 'hist');
hold on;
line([0 0], [0 max(hPFC)*1.2], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
hold off;
xlim(xrange);
ylim([0 max(hPFC)*1.2]);
title('dlFC');
xlabel(sprintf('Fraction change\nafter target'));

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print('12', '-dpdf');

%% ANOVA for post-target change
pt_A1   = mA1_pta-mA1_ptp;
pt_dPEG = mdPEG_pta-mdPEG_ptp;
pt_VPr  = mVPr_pta-mVPr_ptp;
pt_PFC  = mPFC_pta-mPFC_ptp;

pt_mat = NaN(max([length(pt_A1), length(pt_dPEG), length(pt_VPr), length(pt_PFC)]), 4);
pt_mat(1:length(pt_A1),1) = pt_A1;
pt_mat(1:length(pt_dPEG),2) = pt_dPEG;
pt_mat(1:length(pt_VPr),3) = pt_VPr;
pt_mat(1:length(pt_PFC),4) = pt_PFC;

[pt_p, pt_table, pt_stats] = kruskalwallis(pt_mat);
figure;
c = multcompare(pt_stats);%, 'ctype', 'dunn-sidak');


%% figure 13 : change in reference vs change in target scatter plot
if normalize_each == 0
  xrange = [-10 10];
end
if popnorm == 1
  xrange = [-5 5];
else
  xrange = [-20 20];
end


figure;

subplot(2,2,1);
scatter(A1_Dapr,A1_Dapt);
xlim(xrange);ylim(xrange);
line(xrange,[ 0 0]);
line([0 0], xrange);
line(xrange, xrange, 'LineStyle', '--', 'Color', 'k');
title('A1');
xlabel('Change in Reference response');
ylabel('Change in Target response');

subplot(2,2,2);
scatter(dPEG_Dapr,dPEG_Dapt);
xlim(xrange);ylim(xrange);
line(xrange,[ 0 0]);
line([0 0], xrange);
line(xrange, xrange, 'LineStyle', '--', 'Color', 'k');
title('dPEG');
xlabel('Change in Reference response');
ylabel('Change in Target response');

subplot(2,2,3);
scatter(VPr_Dapr,VPr_Dapt);
xlim(xrange);ylim(xrange);
line(xrange,[ 0 0]);
line([0 0], xrange);
line(xrange, xrange, 'LineStyle', '--', 'Color', 'k');
title('VPr');
xlabel('Change in Reference response');
ylabel('Change in Target response');

subplot(2,2,4);
scatter(PFC_Dapr,PFC_Dapt);
xlim(xrange);ylim(xrange);
hold on;
line(xrange,[ 0 0]);
line([0 0], xrange);
line(xrange, xrange, 'LineStyle', '--', 'Color', 'k');
title('PFC');
xlabel('Change in Reference response');
ylabel('Change in Target response');



% %% figure XX : comparison of spontaneous activity
% 
% figure;
% set(gcf, 'Position', [398 242 431 823]);
% subplot(5,2,1);
% hist(A1p_bsl, 20);
% xlim([0 60]);
% title('Passive');
% 
% subplot(5,2,2);
% hist(A1a_bsl, 20);
% xlim([0 60]);
% title('Active');
% 
% subplot(5,2,3);
% hist(dPEGp_bsl, 20);
% xlim([0 60]);
% 
% subplot(5,2,4);
% hist(dPEGa_bsl, 20);
% xlim([0 60]);
% 
% subplot(5,2,5);
% hist(VPrp_bsl, 20);
% xlim([0 60]);
% 
% subplot(5,2,6);
% hist(VPra_bsl, 20);
% xlim([0 60]);
% 
% subplot(5,2,7);
% hist(PFCp_bsl, 20);
% xlim([0 60]);
% 
% subplot(5,2,8);
% hist(PFCa_bsl, 20);
% xlim([0 60]);
% 
% subplot(5,1,5);
% bar([nanmean(A1p_bsl) mean(A1a_bsl); mean(dPEGp_bsl) mean(dPEGa_bsl) ; mean(VPrp_bsl) mean(VPra_bsl) ; nanmean(PFCp_bsl) mean(PFCa_bsl)]);
% set(gca, 'XTickLabel', {'A1', 'dPEG', 'VPr', 'dlFC'});
% 
% suptitle('Spontaneous activity - Click discrimination');
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print('13', '-dpdf');

return;
% keyboard;

%% try to correlate behavioral performance with off response

% defined post reference responses as mAREA_prp and _pra (passive and active)
% post target responses are mAREA_ptp and _pta
% post stim window are ptwin1 and ptwin2 (for CLK from2.1 to 2.8 s)

%post_target_index = (mVPr_pta-mVPr_pra)-(mVPr_ptp-mVPr_prp);
post_target_index = (mVPr_pta-mVPr_ptp); %-(mVPr_ptp-mVPr_prp);

mfiles = cell(1,length(VPr_cells_used));
DR = zeros(1,length(VPr_cells_used));
HR = zeros(1,length(VPr_cells_used));
SR = zeros(1,length(VPr_cells_used));

% need to collect behavioral data for each cellid
for n = 1 : length(VPr_cells_used)
  cellinfo = dbgetscellfile('cellid', VPr_cells_used{n}, 'runclass', runclass);
  
  if strcmp(runclass, 'PTD') && isempty(cellinfo)
    cellinfo = dbgetscellfile('cellid', VPr_cells_used{n}, 'runclass', 'CCH');
  end
    
  for nn = 1 : length(cellinfo)
    if strcmp(cellinfo(nn).behavior, 'active')
      mfiles{n} = [cellinfo(nn).stimpath cellinfo(nn).stimfile, '.m'];
    end
  end
end

for n = 1 : length(VPr_cells_used)
  LoadMFile(mfiles{n});
  DR(n) = exptparams.Performance(end).DiscriminationRate;
  HR(n) = exptparams.Performance(end).HitRate;
  SR(n) = exptparams.Performance(end).SafeRate;
end

keyboard;
%% plot results of post-target response vs behavior performance
behavior_measure = DR;

figure;scatter(post_target_index, behavior_measure)
%xlim([-5 5]);
ylim([0 100]);
line([0 0], [0 100], 'LineStyle', '--', 'Color', 'k');

figure; bar([median(post_target_index(behavior_measure<25)), median(post_target_index(behavior_measure>25 & behavior_measure<50)), median(post_target_index(behavior_measure>50 & behavior_measure<75)), median(post_target_index(behavior_measure>75))])

groups = zeros(length(post_target_index), 1);
groups(behavior_measure<25)          = 1;
groups(behavior_measure>=25 & behavior_measure<50) = 2;
groups(behavior_measure>=50 & behavior_measure<75) = 3;
groups(behavior_measure>=75)         = 4;



figure; boxplot(post_target_index, groups, 'orientation', 'vertical', 'notch', 'on')
grid ON;

%% post-target changes (absolute value)
VPr_apt = abs(post_target_index);
[R, P] = corrcoef(behavior_measure, VPr_apt)

figure; scatter(DR,VPr_apt);
lslh = lsline;
set(lslh, 'Color', 'r');
ms = sprintf('DR vs absolute change in post target firing rate, p = %.3f', P(1,2));
title(ms);
xlabel('DR');
ylabel('abs (target active - target passive)');
text(70,2.8, sprintf('R = %.4f\np = % .4f', R(1,2), P(1,2)));

%% post-target changes (ALL)
figure; scatter(DR,post_target_index);
[R, P] = corrcoef(behavior_measure, post_target_index)
lslh = lsline;
set(lslh, 'Color', 'r');
ms = sprintf('DR vs change in post target firing rate, p = %.3f', P(1,2));
title(ms);
xlabel('DR');
ylabel('target active - target passive');
text(70,2.8, sprintf('R = %.4f\np = % .4f', R(1,2), P(1,2)));

%% post-target enhancements ONLY (positive change)
pos_pti = post_target_index(post_target_index>0);
pos_DR  = DR(post_target_index>0);

figure; scatter(pos_DR,pos_pti);
[R, P] = corrcoef(pos_DR, pos_pti)
lslh = lsline;
set(lslh, 'Color', 'r');
ms = sprintf('DR vs enhancement in post target firing rate, p = %.3f', P(1,2));
title(ms);
xlabel('DR');
ylabel('target active - target passive');
text(70,2.8, sprintf('R = %.4f\np = % .4f', R(1,2), P(1,2)));

%% post-target inhibitions ONLY
neg_pti = post_target_index(post_target_index<=0);
neg_DR  = DR(post_target_index<=0);

figure; scatter(neg_DR,neg_pti);
[R, P] = corrcoef(neg_DR, neg_pti)
lslh = lsline;
set(lslh, 'Color', 'r');
ms = sprintf('DR vs negative change in post target firing rate, p = %.3f', P(1,2));
title(ms);
xlabel('DR');
ylabel('target active - target passive');
text(70,0.5, sprintf('R = %.4f\np = % .4f', R(1,2), P(1,2)));
ylim([-2.5 1]);

%% merge with PTD data
load ('./after_target_resp_behavior/PTD_post_target.mat');
post_target_index = [post_target_index ; PTD_post_target_index];
DR = [DR, PTD_DR];
behavior_measure = [behavior_measure, PTD_DR];
VPr_apt = [VPr_apt; PTD_VPr_apt];

%%
figure; bar([median(post_target_index(behavior_measure<50)),  median(post_target_index(behavior_measure>=50))])

figure; 
chanced= histc(post_target_index(groups==1), [-6:.5:6])
plot([-6:0.5:6], chanced./max(chanced)); hold on


goodd= histc(post_target_index(groups==3), [-6:1:6])
plot([-6:1:6], goodd./max(goodd), 'r');


%% new plots separating UP and DOWN click discrim directions

% UP

% plot 1: passive vs. active

figure;
set(gcf, 'Position', [200 52 509 739]);
if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.4 1.3];
elseif popnorm == 0 && normalize_each == 0
  yrange = [-5 15];
  
end

% A1

subplot(4,2,1);
shadedErrorBar(A1_t, mA1rpU, seA1rpU, '--b', 0);
hold on;
shadedErrorBar(A1_t, mA1raU, seA1raU, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1: Reference pass vs act');

subplot(4,2,2);
shadedErrorBar(A1_t, mA1tpU, seA1tpU, '--r', 0);
hold on;
shadedErrorBar(A1_t, mA1taU, seA1taU, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(A1_cells_usedU)))]);

% dPEG

subplot(4,2,3);
shadedErrorBar(dPEG_t, mdPEGrpU, sedPEGrpU, '--b', 0);
hold on;
shadedErrorBar(dPEG_t, mdPEGraU, sedPEGraU, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG: Reference pass vs act');

subplot(4,2,4);
shadedErrorBar(dPEG_t, mdPEGtpU, sedPEGtpU, '--r', 0);
hold on;
shadedErrorBar(dPEG_t, mdPEGtaU, sedPEGtaU, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(dPEG_cells_usedU)))]);




% VPr
subplot(4,2,5);
shadedErrorBar(VPr_t, mVPrrpU, seVPrrpU, '--b', 0);
hold on;
shadedErrorBar(VPr_t, mVPrraU, seVPrraU, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr: Reference pass vs act');

subplot(4,2,6);
shadedErrorBar(VPr_t, mVPrtpU, seVPrtpU, '--r', 0);
hold on;
shadedErrorBar(VPr_t, mVPrtaU, seVPrtaU, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(VPr_cells_usedU)))]);

% dlFC
subplot(4,2,7);
shadedErrorBar(PFC_t, mPFCrpU, sePFCrpU, '--b', 0);
hold on;
shadedErrorBar(PFC_t, mPFCraU, sePFCraU, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(ylim);
title('dlFC: Reference pass vs act');

subplot(4,2,8);
shadedErrorBar(PFC_t, mPFCtpU, sePFCtpU, '--r', 0);
hold on;
shadedErrorBar(PFC_t, mPFCtaU, sePFCtaU, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dlFC: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(PFC_cells_usedU)))]);

suptitle(sprintf('CRD, passive vs. active\nLow reference vs. high target click rate'));



% plot 2: Reference vs Target

figure;
set(gcf, 'Position', [300 52 509 739]);
set(gcf, 'Color', 'w');
if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.4 1.3];
elseif ~popnorm && ~normalize_each
  yrange = [-20 20];
end

% A1
subplot(4,2,1);
hp1 = shadedErrorBar(A1_t, mA1rpU, seA1rpU, '--b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(A1_t, mA1tpU, seA1tpU, '--r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('\bfPassive', 'FontSize', 14);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');

subplot(4,2,2);
hp1 = shadedErrorBar(A1_t, mA1raU, seA1raU, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(A1_t, mA1taU, seA1taU, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('\bfBehavior', 'FontSize', 14);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(A1_cells_usedU)))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);

% dPEG
if svd_subset == 0
  subplot(4,2,3);
  hp1 = shadedErrorBar(dPEG_t, mdPEGrpU, sedPEGrpU, '--b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(dPEG_t, mdPEGtpU, sedPEGtpU, '--r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  
  subplot(4,2,4);
  hp1 =  shadedErrorBar(dPEG_t, mdPEGraU, sedPEGraU, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  hold on;
  hp1 = shadedErrorBar(dPEG_t, mdPEGtaU, sedPEGtaU, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  title('dPEG: Active ref vs tar');
  text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(dPEG_cells_usedU)))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);  
end

%svd dPEG
if svd_subset == 1
  subplot(4,2,3);
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGrpU, sedPEGrpU, '--b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGtpU, sedPEGtpU, '--r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  
  subplot(4,2,4);
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGraU, sedPEGraU, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  hold on;
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGtaU, sedPEGtaU, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(dPEG_cells_usedU)))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);  
end

% VPr

subplot(4,2,5);
hp1 = shadedErrorBar(VPr_t, mVPrrpU, seVPrrpU, '--b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(VPr_t, mVPrtpU, seVPrtpU, '--r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');

subplot(4,2,6);
hp1 = shadedErrorBar(VPr_t, mVPrraU, seVPrraU, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(VPr_t, mVPrtaU, seVPrtaU, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(VPr_cells_usedU)))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);

% dlFC
subplot(4,2,7);
rh = shadedErrorBar(PFC_t, mPFCrpU, sePFCrpU, '--b', 0);
set(rh.edge, 'LineStyle', 'none');
set(rh.mainLine, 'LineWidth', 2);
hold on;
th = shadedErrorBar(PFC_t, mPFCtpU, sePFCtpU, '--r',   0);
set(rh.edge, 'LineStyle', 'none');
set(th.mainLine, 'LineWidth', 2);
legend([rh.mainLine th.mainLine], 'Reference', 'Target');
legend BOXOFF;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');
xlabel('Time (s)');

subplot(4,2,8);
hp1 = shadedErrorBar(PFC_t, mPFCraU, sePFCraU, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(PFC_t, mPFCtaU, sePFCtaU, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(PFC_cells_usedU)))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);
xlabel('Time (s)');

suptitle(sprintf('CRD, reference vs. target\nLow reference vs. high target click rate'));





%% Plot 5: Naive vs Trained figure
if add_naive == 1
  
  if only_pos_tar == -1
    yrange([-0.5 0.5]);
  else
    yrange =[-0.2 0.6];
  end
  
  if popnorm
    yrange = [-1 1.3];
  elseif ~popnorm && ~normalize_each
    yrange = [-20 20];
  end
  
  figure;
  set(gcf, 'Position', [21 48 693 741]);
  set(gcf, 'Color', 'w');  
  % A1

  subplot(3,3,1);
  h5 = shadedErrorBar(nA1_t, mnA1rpU, senA1rpU, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(nA1_t, mnA1tpU, senA1tpU, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  ylim(yrange);
  xlim([0.9 2.2]);
  hold off;
  title('\bfNaive Passive CRD');
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(nA1_cells_usedU)))]);
  set(gca, 'box', 'off');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfA1', 'FontSize', 12);
  set(gca, 'XTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,2);
  h5 = shadedErrorBar(A1_t, mA1rpU, seA1rpU, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(A1_t, mA1tpU, seA1tpU, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold off;
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  title('\bfTrained Passive CRD');
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(unique(A1_cells_usedU)))]);
  
  subplot(3,3,3);
  h5 = shadedErrorBar(A1_t, mA1raU, seA1raU, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(A1_t, mA1taU, seA1taU, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([1 2.2]);
  ylim(yrange);
  title('\bfBehaving CRD');
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');  
  
  
  % dPEG
  if svd_subset == 1
    t = svd_dPEG_t;
  else
    t = dPEG_t;
  end
  

  subplot(3,3,4);
  h5 = shadedErrorBar(ndPEG_t, mndPEGrpU, sendPEGrpU, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(ndPEG_t, mndPEGtpU, sendPEGtpU, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  ylim(yrange);
  xlim([.9 2.2]);
  hold off;
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(ndPEG_cells_usedU)))]);
  set(gca, 'box', 'off');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfdPEG', 'FontSize', 12);
  set(gca, 'XTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,5);
  h5 = shadedErrorBar(t, mdPEGrpU, sedPEGrpU, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(t, mdPEGtpU, sedPEGtpU, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(unique(dPEG_cells_usedU)))]);
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,6);
  h5 = shadedErrorBar(t, mdPEGraU, sedPEGraU, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  ylim(yrange);
  hold on;
  h5 = shadedErrorBar(t, mdPEGtaU, sedPEGtaU, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);  
  hold off;
  xlim([1 2.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top'); 
  
  
  
  % VPr
  subplot(3,3,7);
  h5 = shadedErrorBar(nVPr_t, mnVPrrpU, senVPrrp, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(nVPr_t, mnVPrtpU, senVPrtp, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  xlim([.9 2.2]);
  ylim(yrange);
  hold off;
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(nVPr_cells_usedU)))]);
  set(gca, 'box', 'off');
  set(gca, 'Layer', 'top');
  xlabel('Time (s)');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfVPr', 'FontSize', 12);
  
  subplot(3,3,8);
  h5 = shadedErrorBar(VPr_t, mVPrrpU, seVPrrpU, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(VPr_t, mVPrtpU, seVPrtpU, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);  
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(unique(VPr_cells_usedU)))]);
  set(gca, 'box', 'off');
  xlabel('Time (s)');
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,9);
  h5 = shadedErrorBar(VPr_t, mVPrraU, seVPrraU, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(VPr_t, mVPrtaU, seVPrtaU, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  set(gca, 'box', 'off');
  xlabel('Time (s)');
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  suptitle(sprintf('CRD, reference vs. target\nLow reference vs. high target click rate'));
  
end







% DOWN

% plot 1: passive vs. active

figure;
set(gcf, 'Position', [200 52 509 739]);
if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.4 1.3];
elseif popnorm == 0 && normalize_each == 0
  yrange = [-5 15];
  
end

% A1

subplot(4,2,1);
shadedErrorBar(A1_t, mA1rpD, seA1rpD, '--b', 0);
hold on;
shadedErrorBar(A1_t, mA1raD, seA1raD, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1: Reference pass vs act');

subplot(4,2,2);
shadedErrorBar(A1_t, mA1tpD, seA1tpD, '--r', 0);
hold on;
shadedErrorBar(A1_t, mA1taD, seA1taD, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('A1: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(A1_cells_usedD)))]);

% dPEG

subplot(4,2,3);
shadedErrorBar(dPEG_t, mdPEGrpD, sedPEGrpD, '--b', 0);
hold on;
shadedErrorBar(dPEG_t, mdPEGraD, sedPEGraD, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG: Reference pass vs act');

subplot(4,2,4);
shadedErrorBar(dPEG_t, mdPEGtpD, sedPEGtpD, '--r', 0);
hold on;
shadedErrorBar(dPEG_t, mdPEGtaD, sedPEGtaD, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dPEG: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(dPEG_cells_usedD)))]);




% VPr
subplot(4,2,5);
shadedErrorBar(VPr_t, mVPrrpD, seVPrrpD, '--b', 0);
hold on;
shadedErrorBar(VPr_t, mVPrraD, seVPrraD, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr: Reference pass vs act');

subplot(4,2,6);
shadedErrorBar(VPr_t, mVPrtpD, seVPrtpD, '--r', 0);
hold on;
shadedErrorBar(VPr_t, mVPrtaD, seVPrtaD, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('VPr: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(VPr_cells_usedD)))]);

% dlFC
subplot(4,2,7);
%shadedErrorBar(PFC_t, mPFCrpD, sePFCrpD, '--b', 0);
hold on;
%shadedErrorBar(PFC_t, mPFCraD, sePFCraD, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(ylim);
title('dlFC: Reference pass vs act');

subplot(4,2,8);
%shadedErrorBar(PFC_t, mPFCtpD, sePFCtpD, '--r', 0);
hold on;
%shadedErrorBar(PFC_t, mPFCtaD, sePFCtaD, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('dlFC: Target pass vs act');
text(0.1, yrange(2)*.8, ['N=' num2str(length(unique(PFC_cells_usedD)))]);

suptitle(sprintf('CRD, passive vs. active\nHigh reference vs. low target click rate'));



% plot 2: Reference vs Target

figure;
set(gcf, 'Position', [300 52 509 739]);
set(gcf, 'Color', 'w');
if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.4 1.3];
elseif ~popnorm && ~normalize_each
  yrange = [-20 20];
end

% A1
subplot(4,2,1);
hp1 = shadedErrorBar(A1_t, mA1rpD, seA1rpD, '--b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(A1_t, mA1tpD, seA1tpD, '--r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('\bfPassive', 'FontSize', 14);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');

subplot(4,2,2);
hp1 = shadedErrorBar(A1_t, mA1raD, seA1raD, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(A1_t, mA1taD, seA1taD, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
title('\bfBehavior', 'FontSize', 14);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(A1_cells_usedD)))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);

% dPEG
if svd_subset == 0
  subplot(4,2,3);
  hp1 = shadedErrorBar(dPEG_t, mdPEGrpD, sedPEGrpD, '--b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(dPEG_t, mdPEGtpD, sedPEGtpD, '--r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  
  subplot(4,2,4);
  hp1 =  shadedErrorBar(dPEG_t, mdPEGraD, sedPEGraD, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  hold on;
  hp1 = shadedErrorBar(dPEG_t, mdPEGtaD, sedPEGtaD, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  title('dPEG: Active ref vs tar');
  text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(dPEG_cells_usedD)))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);  
end

%svd dPEG
if svd_subset == 1
  subplot(4,2,3);
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGrpD, sedPEGrpD, '--b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGtpD, sedPEGtpD, '--r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  
  subplot(4,2,4);
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGraD, sedPEGraD, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  hold on;
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGtaD, sedPEGtaD, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);  
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.4 3.2]);
  ylim(yrange);
  text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(dPEG_cells_usedD)))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);  
end

% VPr

subplot(4,2,5);
hp1 = shadedErrorBar(VPr_t, mVPrrpD, seVPrrpD, '--b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(VPr_t, mVPrtpD, seVPrtpD, '--r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');

subplot(4,2,6);
hp1 = shadedErrorBar(VPr_t, mVPrraD, seVPrraD, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(VPr_t, mVPrtaD, seVPrtaD, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(VPr_cells_usedD)))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);

% dlFC
subplot(4,2,7);
%rh = shadedErrorBar(PFC_t, mPFCrpD, sePFCrpD, '--b', 0);
set(rh.edge, 'LineStyle', 'none');
set(rh.mainLine, 'LineWidth', 2);
hold on;
%th = shadedErrorBar(PFC_t, mPFCtpD, sePFCtpD, '--r',   0);
set(rh.edge, 'LineStyle', 'none');
set(th.mainLine, 'LineWidth', 2);
legend([rh.mainLine th.mainLine], 'Reference', 'Target');
legend BOXOFF;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');
xlabel('Time (s)');

subplot(4,2,8);
%hp1 = shadedErrorBar(PFC_t, mPFCraD, sePFCraD, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
%hp1 = shadedErrorBar(PFC_t, mPFCtaD, sePFCtaD, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.4 3.2]);
ylim(yrange);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(unique(PFC_cells_usedD)))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);
xlabel('Time (s)');

suptitle(sprintf('CRD, reference vs. target\nHigh reference vs. low target click rate'));





%% Plot 5: Naive vs Trained figure
if add_naive == 1
  
  if only_pos_tar == -1
    yrange([-0.5 0.5]);
  else
    yrange =[-0.2 0.6];
  end
  
  if popnorm
    yrange = [-1 1.3];
  elseif ~popnorm && ~normalize_each
    yrange = [-20 20];
  end
  
  figure;
  set(gcf, 'Position', [21 48 693 741]);
  set(gcf, 'Color', 'w');  
  % A1

  subplot(3,3,1);
  h5 = shadedErrorBar(nA1_t, mnA1rpD, senA1rpD, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(nA1_t, mnA1tpD, senA1tpD, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  ylim(yrange);
  xlim([0.9 2.2]);
  hold off;
  title('\bfNaive Passive CRD');
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(nA1_cells_usedD)))]);
  set(gca, 'box', 'off');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfA1', 'FontSize', 12);
  set(gca, 'XTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,2);
  h5 = shadedErrorBar(A1_t, mA1rpD, seA1rpD, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(A1_t, mA1tpD, seA1tpD, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold off;
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  title('\bfTrained Passive CRD');
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(unique(A1_cells_usedD)))]);
  
  subplot(3,3,3);
  h5 = shadedErrorBar(A1_t, mA1raD, seA1raD, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(A1_t, mA1taD, seA1taD, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([1 2.2]);
  ylim(yrange);
  title('\bfBehaving CRD');
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');  
  
  
  % dPEG
  if svd_subset == 1
    t = svd_dPEG_t;
  else
    t = dPEG_t;
  end
  

  subplot(3,3,4);
  h5 = shadedErrorBar(ndPEG_t, mndPEGrpD, sendPEGrpD, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(ndPEG_t, mndPEGtpD, sendPEGtpD, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  ylim(yrange);
  xlim([.9 2.2]);
  hold off;
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(ndPEG_cells_usedD)))]);
  set(gca, 'box', 'off');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfdPEG', 'FontSize', 12);
  set(gca, 'XTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,5);
  h5 = shadedErrorBar(t, mdPEGrpD, sedPEGrpD, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(t, mdPEGtpD, sedPEGtpD, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(unique(dPEG_cells_usedD)))]);
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,6);
  h5 = shadedErrorBar(t, mdPEGraD, sedPEGraD, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  ylim(yrange);
  hold on;
  h5 = shadedErrorBar(t, mdPEGtaD, sedPEGtaD, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);  
  hold off;
  xlim([1 2.2]);
  ylim(yrange);
  set(gca, 'box', 'off');
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top'); 
  
  
  
  % VPr
  subplot(3,3,7);
  h5 = shadedErrorBar(nVPr_t, mnVPrrpD, senVPrrpD, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(nVPr_t, mnVPrtpD, senVPrtpD, 'r', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  xlim([.9 2.2]);
  ylim(yrange);
  hold off;
  text(1.8, yrange(2)*.8, ['\bfN=' num2str(length(unique(nVPr_cells_usedD)))]);
  set(gca, 'box', 'off');
  set(gca, 'Layer', 'top');
  xlabel('Time (s)');
  ylabel('Normalized firing rate');
  text(0.2, 0.4, '\bfVPr', 'FontSize', 12);
  
  subplot(3,3,8);
  h5 = shadedErrorBar(VPr_t, mVPrrpD, seVPrrpD, 'b',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(VPr_t, mVPrtpD, seVPrtpD, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);  
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  text(1.9, yrange(2)*.8, ['\bfN=' num2str(length(unique(VPr_cells_usedD)))]);
  set(gca, 'box', 'off');
  xlabel('Time (s)');
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  subplot(3,3,9);
  h5 = shadedErrorBar(VPr_t, mVPrraD, seVPrraD, 'b', 0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  hold on;
  h5 = shadedErrorBar(VPr_t, mVPrtaD, seVPrtaD, 'r',   0);
  set(h5.edge, 'LineStyle', 'none');
  set(h5.mainLine, 'LineWidth', 2);  
  line([1.25 1.25],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  ylim(yrange);
  xlim([1 2.2]);
  set(gca, 'box', 'off');
  xlabel('Time (s)');
  set(gca, 'YTickLabel', []);
  set(gca, 'Layer', 'top');
  
  suptitle(sprintf('CRD, reference vs. target\nHigh reference vs. low target click rate'));
  
end

return;



%% organize and save variables for PTD - CLK comparisons

CLK_A1_cells_used = A1_cells_used;
CLK_A1_TEI = A1_Dtra - A1_Dtrp;
CLK_A1_Dtrp = A1_Dtrp;
CLK_A1_Dtra = A1_Dtra;
CLK_A1_at = A1_at;
CLK_A1_Dapr = A1_Dapr;
CLK_A1_Dapt = A1_Dapt;

CLK_dPEG_cells_used = dPEG_cells_used;
CLK_dPEG_TEI = dPEG_Dtra - dPEG_Dtrp;
CLK_dPEG_Dtrp = dPEG_Dtrp;
CLK_dPEG_Dtra = dPEG_Dtra;
CLK_dPEG_at = dPEG_at;
CLK_dPEG_Dapr = dPEG_Dapr;
CLK_dPEG_Dapt = dPEG_Dapt;

CLK_VPr_cells_used = VPr_cells_used;
CLK_VPr_TEI = VPr_Dtra - VPr_Dtrp;
CLK_VPr_Dtrp = VPr_Dtrp;
CLK_VPr_Dtra = VPr_Dtra;
CLK_VPr_at = VPr_at;
CLK_VPr_Dapr = VPr_Dapr;
CLK_VPr_Dapt = VPr_Dapt;

CLK_PFC_cells_used = PFC_cells_used;
CLK_PFC_TEI = PFC_Dtra - PFC_Dtrp;
CLK_PFC_Dtrp = PFC_Dtrp;
CLK_PFC_Dtra = PFC_Dtra;
CLK_PFC_at = PFC_at;
CLK_PFC_Dapr = PFC_Dapr;
CLK_PFC_Dapt = PFC_Dapt;

save('CLK_scatvals.mat', 'CLK_*');


%% TRY to come up with an index of change
% Jonathan's idea of calculating the ratio if the change in the area under the curves
% of figure 4 in the ms

% during 0.1 to 0.45 from onset (scatwin)
anawin  = [51:96]; %from 1.25s to 2.4s (click part)
A1enh = (trapz(A1_t(anawin), mA1ta(anawin)) - trapz(A1_t(anawin), mA1ra(anawin))) / (trapz(A1_t(anawin), mA1tp(anawin)) - trapz(A1_t(anawin), mA1rp(anawin)));
dPEGenh = (trapz(dPEG_t(anawin), mdPEGta(anawin)) - trapz(dPEG_t(anawin), mdPEGra(anawin))) / (trapz(dPEG_t(anawin), mdPEGtp(anawin)) - trapz(dPEG_t(anawin), mdPEGrp(anawin)));
VPrenh = (trapz(VPr_t(anawin), mVPrta(anawin)) - trapz(VPr_t(anawin), mVPrra(anawin))) / (trapz(VPr_t(anawin), mVPrtp(anawin)) - trapz(VPr_t(anawin), mVPrrp(anawin)));
PFCenh = (trapz(PFC_t(anawin), mPFCta(anawin)) - trapz(PFC_t(anawin), mPFCra(anawin))) / abs((trapz(PFC_t(anawin), mPFCtp(anawin)) - trapz(PFC_t(anawin), mPFCrp(anawin))));
figure;
plot([A1enh, dPEGenh, VPrenh, PFCenh])
hold on;
plot([A1enh, dPEGenh, VPrenh, PFCenh], 'ro');


% from onset to offset+.4s, separating passive and active
anawin  = [51:85]; %from 1.25s to 2.4s (click part)
% lastt   = A1_t(anawin(end));

A1enhP = trapz(A1_t(anawin), mA1tp(anawin)) - trapz(A1_t(anawin), mA1rp(anawin));
A1enhA = trapz(A1_t(anawin), mA1ta(anawin)) - trapz(A1_t(anawin), mA1ra(anawin));
dPEGenhP = trapz(dPEG_t(anawin), mdPEGtp(anawin)) - trapz(dPEG_t(anawin), mdPEGrp(anawin));
dPEGenhA = trapz(dPEG_t(anawin), mdPEGta(anawin)) - trapz(dPEG_t(anawin), mdPEGra(anawin));
VPrenhP = trapz(VPr_t(anawin), mVPrtp(anawin)) - trapz(VPr_t(anawin), mVPrrp(anawin));
VPrenhA = trapz(VPr_t(anawin), mVPrta(anawin)) - trapz(VPr_t(anawin), mVPrra(anawin));
PFCenhP = trapz(PFC_t(anawin), mPFCtp(anawin)) - trapz(PFC_t(anawin), mPFCrp(anawin));
PFCenhA = trapz(PFC_t(anawin), mPFCta(anawin)) - trapz(PFC_t(anawin), mPFCra(anawin));

% A1enhP   = A1enhP   / lastt;
% A1enhA   = A1enhA   / lastt;
% dPEGenhP = dPEGenhP / lastt;
% dPEGenhA = dPEGenhA / lastt;
% VPrenhP  = VPrenhP  / lastt;
% VPrenhA  = VPrenhA  / lastt;
% PFCenhP  = PFCenhP  / lastt;
% PFCenhA  = PFCenhA  / lastt;

figure;
subplot(2,1,1);
plot([A1enhP, dPEGenhP, VPrenhP, PFCenhP], 'b')
hold on;
plot([A1enhP, dPEGenhP, VPrenhP, PFCenhP], 'bo');
plot([A1enhA, dPEGenhA, VPrenhA, PFCenhA], 'r')
hold on;
plot([A1enhA, dPEGenhA, VPrenhA, PFCenhA], 'ro');

subplot(2,1,2);
bar([A1enhA/A1enhP, dPEGenhA/dPEGenhP, VPrenhA/VPrenhP, PFCenhA/PFCenhP]);

% from onset to offset+.4s, ref passive vs target active (BAD IDEA)
anawin  = [51:85]; %from 1.25s to 2.4s (click part)

lastt   = A1_t(anawin(end));

A1enh = trapz(A1_t(anawin), mA1ta(anawin)) - trapz(A1_t(anawin), mA1rp(anawin));
dPEGenh = trapz(dPEG_t(anawin), mdPEGta(anawin)) - trapz(dPEG_t(anawin), mdPEGrp(anawin));
VPrenh = trapz(VPr_t(anawin), mVPrta(anawin)) - trapz(VPr_t(anawin), mVPrrp(anawin));
PFCenh = trapz(PFC_t(anawin), mPFCta(anawin)) - trapz(PFC_t(anawin), mPFCrp(anawin));

A1enh   = A1enh   / lastt;
dPEGenh = dPEGenh / lastt;
VPrenh  = VPrenh  / lastt;
PFCenh  = PFCenh  / lastt;

figure;
bar([A1enh, dPEGenh, VPrenh, PFCenh]);
