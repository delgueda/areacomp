% script to play with normalizations, statistics, and plotting of my data
% Tone detection
% DEG 6-May-2017

% using data from SC_get_meanPSTH.m and
% naive_MU_context_analysis.m

load('./30Hz-nomodfilt/PTD_PSTHs.mat');
%load('./with_post_30Hz/PTD_PSTHs.mat');
load('./30Hz-nomodfilt/naive_TONE_PSTH.mat'); % SU
%load('./30Hz/naive_TONE_PSTH.mat'); % MUs

latstat = 0;

if 0 % correct number of cells to the list we are using in the paper
  load ('./with_post_30Hz/official_cell_list.mat');
  A1_keep = ismember(A1_cells_used, neurons.A1);
  A1_cells_used = A1_cells_used(A1_keep);
  A1_FTC        = A1_FTC(A1_keep,:);
  A1_Ractive    = A1_Ractive(A1_keep,:);
  A1_Rpassive   = A1_Rpassive(A1_keep,:);
  A1_Rpost1     = A1_Rpost1(A1_keep,:);
  A1_Rpost2     = A1_Rpost2(A1_keep,:);
  A1_TOR        = A1_TOR(A1_keep,:);
  A1_Tactive    = A1_Tactive(A1_keep,:);
  A1_Tpassive   = A1_Tpassive(A1_keep,:);
  A1_Tpost1     = A1_Tpost1(A1_keep,:);
  A1_Tpost2     = A1_Tpost2(A1_keep,:);
  A1_ftc_mods   = A1_ftc_mods(A1_keep);
  A1_params     = A1_params(A1_keep);
  A1_sig_eff    = A1_sig_eff(A1_keep);
  A1_tor_mods   = A1_tor_mods(A1_keep);
  
  dPEG_keep = ismember(dPEG_cells_used, neurons.dPEG);
  dPEG_cells_used = dPEG_cells_used(dPEG_keep);
  dPEG_FTC        = dPEG_FTC(dPEG_keep,:);
  dPEG_Ractive    = dPEG_Ractive(dPEG_keep,:);
  dPEG_Rpassive   = dPEG_Rpassive(dPEG_keep,:);
  dPEG_Rpost1     = dPEG_Rpost1(dPEG_keep,:);
  dPEG_Rpost2     = dPEG_Rpost2(dPEG_keep,:);
  dPEG_TOR        = dPEG_TOR(dPEG_keep,:);
  dPEG_Tactive    = dPEG_Tactive(dPEG_keep,:);
  dPEG_Tpassive   = dPEG_Tpassive(dPEG_keep,:);
  dPEG_Tpost1     = dPEG_Tpost1(dPEG_keep,:);
  dPEG_Tpost2     = dPEG_Tpost2(dPEG_keep,:);
  dPEG_ftc_mods   = dPEG_ftc_mods(dPEG_keep);
  dPEG_params     = dPEG_params(dPEG_keep);
  dPEG_sig_eff    = dPEG_sig_eff(dPEG_keep);
  dPEG_tor_mods   = dPEG_tor_mods(dPEG_keep);
end

% normalization settings
subtract_baseline = 1; % PAPER DEFAULT 1
normalize_each    = 0; % normalize by each cell
   indep_norm     = 0; % 1:independently normalize passives and actives; 0: normalize pass and act by same value (max across all)
popnorm           = 1; % PAPER DEFAULT 1

% selection rules
limit_single_tar  = 0; %1: only 1 tar, 0:any -1:exclude 1 tar
only_pos_tar      = 0; %1:only +, 0: all, -1:exclude + resp
   pos_tar_act       = 1; %1:select positive targets in active, 0: passive target
svd_subset        = 0;

plot_ctx_figures = 0; % 1 for plotting the task-context figures
scatwin = [10:20]; % bins to use for calculating the target-reference contrast in scatter and hist figures
% 28nov2017 : default 10:20

A1p_bsl = [];
A1a_bsl = [];
dPEGp_bsl = [];
dPEGa_bsl = [];
VPrp_bsl = [];
VPra_bsl = [];
PFCp_bsl = [];
PFCa_bsl = [];

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
    
    bsl = A1_TOR(i,1:tonset);
    A1_TOR(i,:) = A1_TOR(i,:) - nanmean(bsl);
    
    bsl = A1_FTC(i,1:tonset);
    A1_FTC(i,:) = A1_FTC(i,:) - nanmean(bsl);
  end
  
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    
    if indep_norm == 1
      fmaxP = max([A1_Rpassive(i,:) A1_Tpassive(i,:)]);
      fmaxA = max([A1_Ractive(i,:) A1_Tactive(i,:)]);
    elseif indep_norm == 0
      fmaxP = max([A1_Rpassive(i,:) A1_Tpassive(i,:) A1_Ractive(i,:) A1_Tactive(i,:)]);
      fmaxA = fmaxP; 
    end
    
    A1_Rpassive(i,:) = A1_Rpassive(i,:) ./ fmaxP;
    A1_Tpassive(i,:) = A1_Tpassive(i,:) ./ fmaxP;
    A1_Ractive(i,:) = A1_Ractive(i,:)   ./ fmaxA;
    A1_Tactive(i,:) = A1_Tactive(i,:)   ./ fmaxA;
    
    A1_TOR(i,:)     = A1_TOR(i,:)       ./ fmaxP;
    A1_FTC(i,:)     = A1_FTC(i,:)       ./ fmaxP;
  end
  
  
end

if popnorm
  [A1_Rpassive, A1_Ractive, A1_Tpassive, A1_Tactive,A1_FTC, A1_TOR] = popmaxnorm...
    (A1_Rpassive, A1_Ractive, A1_Tpassive, A1_Tactive,A1_FTC, A1_TOR);

end

% naive (tango) data
for i=1:size(nA1_Rpassive, 1)
  % baseline mean removal
  
  tonset = 0.2 * nA1_sr; % presilence was fixed at 0.2s
  %bsl = zeros(1,nbins).*NaN;
  
  if subtract_baseline
    bsl = nA1_Rpassive(i,1:tonset);
    nA1_Rpassive(i,:) = nA1_Rpassive(i,:) - nanmean(bsl);
    
    bsl = nA1_TOR(i,1:tonset);
    nA1_TOR(i,:) = nA1_TOR(i,:) - nanmean(bsl);
    
    bsl = nA1_Tpassive(i,1:tonset);
    nA1_Tpassive(i,:) = nA1_Tpassive(i,:) - nanmean(bsl);
    
    bsl = nA1_FTC(i,1:tonset);
    nA1_FTC(i,:) = nA1_FTC(i,:) - nanmean(bsl);
    
  end
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    fmax = max([nA1_Rpassive(i,:) nA1_Tpassive(i,:) nA1_FTC(i,:) nA1_TOR(i,:)]);
    nA1_Rpassive(i,:) = nA1_Rpassive(i,:) ./ fmax;
    nA1_Tpassive(i,:) = nA1_Tpassive(i,:) ./ fmax;
    
    nA1_TOR(i,:)     = nA1_TOR(i,:)       ./ fmax;
    nA1_FTC(i,:)     = nA1_FTC(i,:)       ./ fmax;
  end
  
end
if popnorm
  [nA1_Rpassive, nA1_Tpassive] = popmaxnorm(nA1_Rpassive, nA1_Tpassive); % should I add FTC and TOR?
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
    
    bsl = dPEG_TOR(i,1:tonset);
    dPEG_TOR(i,:) = dPEG_TOR(i,:) - nanmean(bsl);
    
    bsl = dPEG_FTC(i,1:tonset);
    dPEG_FTC(i,:) = dPEG_FTC(i,:) - nanmean(bsl);
  end
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 1
      fmaxP = max([dPEG_Rpassive(i,:) dPEG_Tpassive(i,:)]);
      fmaxA = max([dPEG_Ractive(i,:) dPEG_Tactive(i,:)]);
    elseif indep_norm == 0
      fmaxP = max([dPEG_Rpassive(i,:) dPEG_Tpassive(i,:) dPEG_Ractive(i,:) dPEG_Tactive(i,:)]);
      fmaxA = fmaxP;
    end
    
    dPEG_Rpassive(i,:) = dPEG_Rpassive(i,:) ./ fmaxP;
    dPEG_Tpassive(i,:) = dPEG_Tpassive(i,:) ./ fmaxP;
    dPEG_Ractive(i,:) = dPEG_Ractive(i,:)   ./ fmaxA;
    dPEG_Tactive(i,:) = dPEG_Tactive(i,:)   ./ fmaxA;
    
    dPEG_TOR(i,:)     = dPEG_TOR(i,:)       ./ fmaxP;
    dPEG_FTC(i,:)     = dPEG_FTC(i,:)       ./ fmaxP;
  end
  
  
end

if popnorm
  [dPEG_Rpassive, dPEG_Ractive, dPEG_Tpassive, dPEG_Tactive,dPEG_FTC, dPEG_TOR] = popmaxnorm...
    (dPEG_Rpassive, dPEG_Ractive, dPEG_Tpassive, dPEG_Tactive,dPEG_FTC, dPEG_TOR);
end

% naive (tango) data
for i=1:size(ndPEG_Rpassive, 1)
  % baseline mean removal
  
  tonset = 0.2 * ndPEG_sr; % presilence was fixed at 0.2s
  %bsl = zeros(1,nbins).*NaN;
  
  if subtract_baseline
    bsl = ndPEG_Rpassive(i,1:tonset);
    ndPEG_Rpassive(i,:) = ndPEG_Rpassive(i,:) - nanmean(bsl);
    
    bsl = ndPEG_TOR(i,1:tonset);
    ndPEG_TOR(i,:) = ndPEG_TOR(i,:) - nanmean(bsl);
    
    bsl = ndPEG_Tpassive(i,1:tonset);
    ndPEG_Tpassive(i,:) = ndPEG_Tpassive(i,:) - nanmean(bsl);
    
    bsl = ndPEG_FTC(i,1:tonset);
    ndPEG_FTC(i,:) = ndPEG_FTC(i,:) - nanmean(bsl);
    
  end
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    fmax = max([ndPEG_Rpassive(i,:) ndPEG_Tpassive(i,:) ndPEG_FTC(i,:) ndPEG_TOR(i,:)]);
    ndPEG_Rpassive(i,:) = ndPEG_Rpassive(i,:) ./ fmax;
    ndPEG_Tpassive(i,:) = ndPEG_Tpassive(i,:) ./ fmax;
    
    ndPEG_TOR(i,:)     = ndPEG_TOR(i,:)       ./ fmax;
    ndPEG_FTC(i,:)     = ndPEG_FTC(i,:)       ./ fmax;
  end
  
end

if popnorm
  [ndPEG_Rpassive, ndPEG_Tpassive] = popmaxnorm(ndPEG_Rpassive, ndPEG_Tpassive);
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
    
    bsl = VPr_TOR(i,1:tonset);
    VPr_TOR(i,:) = VPr_TOR(i,:) - nanmean(bsl);
    
    bsl = VPr_FTC(i,1:tonset);
    VPr_FTC(i,:) = VPr_FTC(i,:) - nanmean(bsl);
  end
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 1
      fmaxP = max([VPr_Rpassive(i,:) VPr_Tpassive(i,:)]);
      fmaxA = max([VPr_Ractive(i,:) VPr_Tactive(i,:)]);
    elseif indep_norm == 0
      fmaxP = max([VPr_Rpassive(i,:) VPr_Tpassive(i,:) VPr_Ractive(i,:) VPr_Tactive(i,:)]);
      fmaxA = fmaxP;
    end
    
    VPr_Rpassive(i,:) = VPr_Rpassive(i,:) ./ fmaxP;
    VPr_Tpassive(i,:) = VPr_Tpassive(i,:) ./ fmaxP;
    VPr_Ractive(i,:) = VPr_Ractive(i,:)   ./ fmaxA;
    VPr_Tactive(i,:) = VPr_Tactive(i,:)   ./ fmaxA;
    
    VPr_TOR(i,:)     = VPr_TOR(i,:)       ./ fmaxP;
    VPr_FTC(i,:)     = VPr_FTC(i,:)       ./ fmaxP;
    
  end
  
  
end

if popnorm
  [VPr_Rpassive, VPr_Ractive, VPr_Tpassive, VPr_Tactive,VPr_FTC, VPr_TOR] = popmaxnorm...
    (VPr_Rpassive, VPr_Ractive, VPr_Tpassive, VPr_Tactive,VPr_FTC, VPr_TOR);
end


% naive (tango) data
for i=1:size(nVPr_Rpassive, 1)
  % baseline mean removal
  
  tonset = 0.2 * nVPr_sr; % presilence was fixed at 0.2s
  %bsl = zeros(1,nbins).*NaN;
  
  if subtract_baseline
    bsl = nVPr_Rpassive(i,1:tonset);
    nVPr_Rpassive(i,:) = nVPr_Rpassive(i,:) - nanmean(bsl);
    
    bsl = nVPr_TOR(i,1:tonset);
    nVPr_TOR(i,:) = nVPr_TOR(i,:) - nanmean(bsl);
    
    bsl = nVPr_Tpassive(i,1:tonset);
    nVPr_Tpassive(i,:) = nVPr_Tpassive(i,:) - nanmean(bsl);
    
    bsl = nVPr_FTC(i,1:tonset);
    nVPr_FTC(i,:) = nVPr_FTC(i,:) - nanmean(bsl);
    
  end
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    fmax = max([nVPr_Rpassive(i,:) nVPr_Tpassive(i,:) nVPr_FTC(i,:) nVPr_TOR(i,:)]);
    nVPr_Rpassive(i,:) = nVPr_Rpassive(i,:) ./ fmax;
    nVPr_Tpassive(i,:) = nVPr_Tpassive(i,:) ./ fmax;
    
    nVPr_TOR(i,:)     = nVPr_TOR(i,:)       ./ fmax;
    nVPr_FTC(i,:)     = nVPr_FTC(i,:)       ./ fmax;
  end
  
end

if popnorm
  [nVPr_Rpassive, nVPr_Tpassive] = popmaxnorm(nVPr_Rpassive, nVPr_Tpassive);
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
    
    bsl = PFC_TOR(i,1:tonset);
    PFC_TOR(i,:) = PFC_TOR(i,:) - nanmean(bsl);
    
    bsl = PFC_FTC(i,1:tonset);
    PFC_FTC(i,:) = PFC_FTC(i,:) - nanmean(bsl);
  end
  
  
  if normalize_each
    % max will be taken from all conditions and sound to each cell
    if indep_norm == 1
      fmaxP = max([PFC_Rpassive(i,:) PFC_Tpassive(i,:)]);
      fmaxA = max([PFC_Ractive(i,:) PFC_Tactive(i,:)]);
    elseif indep_norm == 0
      fmaxP = max([PFC_Rpassive(i,:) PFC_Tpassive(i,:) PFC_Ractive(i,:) PFC_Tactive(i,:)]);
      fmaxA = fmaxP;
    end
    
    PFC_Rpassive(i,:) = PFC_Rpassive(i,:) ./ fmaxP;
    PFC_Tpassive(i,:) = PFC_Tpassive(i,:) ./ fmaxP;
    PFC_Ractive(i,:) = PFC_Ractive(i,:)   ./ fmaxA;
    PFC_Tactive(i,:) = PFC_Tactive(i,:)   ./ fmaxA;
    
    PFC_TOR(i,:)     = PFC_TOR(i,:)       ./ fmaxP;
    PFC_FTC(i,:)     = PFC_FTC(i,:)       ./ fmaxP;
  end
  
end

if popnorm
  [PFC_Rpassive, PFC_Ractive, PFC_Tpassive, PFC_Tactive, PFC_FTC, PFC_TOR] = popmaxnorm...
    (PFC_Rpassive, PFC_Ractive, PFC_Tpassive, PFC_Tactive, PFC_FTC, PFC_TOR);
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
      
      bsl = svd_dPEG_TOR(i,1:tonset);
      svd_dPEG_TOR(i,:) = svd_dPEG_TOR(i,:) - nanmean(bsl);
      
      bsl = svd_dPEG_FTC(i,1:tonset);
      svd_dPEG_FTC(i,:) = svd_dPEG_FTC(i,:) - nanmean(bsl);
    end
    
    if normalize_each
      % max will be taken from all conditions and sound to each cell
      if indep_norm == 1
        fmaxP = max([svd_dPEG_Rpassive(i,:) svd_dPEG_Tpassive(i,:)]);
        fmaxA = max([svd_dPEG_Ractive(i,:) svd_dPEG_Tactive(i,:)]);
      elseif indep_norm == 0
        fmaxP = max([svd_dPEG_Rpassive(i,:) svd_dPEG_Tpassive(i,:) svd_dPEG_Ractive(i,:) svd_dPEG_Tactive(i,:)]);
        fmaxA = fmaxP;
      end
      
      svd_dPEG_Rpassive(i,:) = svd_dPEG_Rpassive(i,:) ./ fmaxP;
      svd_dPEG_Tpassive(i,:) = svd_dPEG_Tpassive(i,:) ./ fmaxP;
      svd_dPEG_Ractive(i,:) = svd_dPEG_Ractive(i,:)   ./ fmaxA;
      svd_dPEG_Tactive(i,:) = svd_dPEG_Tactive(i,:)   ./ fmaxA;
      
      svd_dPEG_TOR(i,:)     = svd_dPEG_TOR(i,:)       ./ fmaxP;
      svd_dPEG_FTC(i,:)     = svd_dPEG_FTC(i,:)       ./ fmaxP;
      
    end
    
    
  end
  
  if popnorm
    [svd_dPEG_Rpassive, svd_dPEG_Ractive, svd_dPEG_Tpassive, svd_dPEG_Tactive, svd_dPEG_FTC, svd_dPEG_TOR] = popmaxnorm...
      (svd_dPEG_Rpassive, svd_dPEG_Ractive, svd_dPEG_Tpassive, svd_dPEG_Tactive, svd_dPEG_FTC, svd_dPEG_TOR);
  end
  
  svd_dPEG_t(end) = []; % there was an extra bin at the end, check why
  
end

%% Selection rules

%%---A1
% selection of cells based on target positive *onset* responses (tpr) during behavior
if only_pos_tar == 1
  if pos_tar_act == 1
    A1tpr = [find(mean(A1_Tactive(:,7:15),2)>0)];
  elseif pos_tar_act == 0
    A1tpr = [find(mean(A1_Tpassive(:,7:15),2)>0)];
  end
  
elseif only_pos_tar ==-1
  if pos_tar_act == 1
    A1tpr = [find(mean(A1_Tactive(:,7:15),2)<=0)];
  elseif pos_tar_act == 0
    A1tpr = [find(mean(A1_Tpassive(:,7:15),2)<=0)];
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
      dPEGtpr = [find(mean(dPEG_Tpassive(:,7:15),2)>0)];
    elseif pos_tar_act
      dPEGtpr = [find(mean(dPEG_Tactive(:,7:15),2)>0)];
    end
    
    %dPEGtpr = find(mean([dPEG_Tactive(:,7:15) dPEG_Tpassive(:,7:15)],2)>0);
  elseif only_pos_tar == -1
    if pos_tar_act == 1
      dPEGtpr = [find(mean(dPEG_Tactive(:,7:15),2)<=0)];
    elseif pos_tar_act == 0
      dPEGtpr = [find(mean(dPEG_Tpassive(:,7:15),2)<=0)];
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
    if pos_tar_act == 1
      dPEGtpr = [find(mean(svd_dPEG_Tactive(:,7:15),2)>0)];
    elseif pos_tar_act == 0
      dPEGtpr = [find(mean(svd_dPEG_Tpassive(:,7:15),2)>0)];
    end
    %dPEGtpr = find(mean([svd_dPEG_Tactive(:,7:15) svd_dPEG_Tpassive(:,7:15)],2)>0);
  elseif only_pos_tar == -1
    if pos_tar_act == 0
      dPEGtpr = [find(mean(svd_dPEG_Tpassive(:,7:15),2)<=0)];
    elseif pos_tar_act == 1
      dPEGtpr = [find(mean(svd_dPEG_Tactive(:,7:15),2)<=0)];
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
  if pos_tar_act == 0
    VPrtpr = [find(mean(VPr_Tpassive(:,7:15),2)>0)];
  elseif pos_tar_act == 1
    VPrtpr = [find(mean(VPr_Tactive(:,7:15),2)>0)];
  end
elseif only_pos_tar == -1
  if pos_tar_act == 0
    VPrtpr = [find(mean(VPr_Tpassive(:,7:15),2)<=0)];
  elseif pos_tar_act == 1
    VPrtpr = [find(mean(VPr_Tactive(:,7:15),2)<=0)];
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
    PFCtpr = [find(mean(PFC_Tactive(:,7:15),2)>0)];
  elseif pos_tar_act == 0
    PFCtpr = [find(mean(PFC_Tpassive(:,7:15),2)>0)];
  end
elseif only_pos_tar == -1
  if pos_tar_act == 0
    PFCtpr = [find(mean(PFC_Tpassive(:,7:15),2)<=0)];
  elseif pos_tar_act == 1
    PFCtpr = [find(mean(PFC_Tactive(:,7:15),2)<=0)];
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
elseif ~popnorm && ~normalize_each
  yrange = [-15 15];
end

% A1
[mA1rp seA1rp] = jackmeanerr(A1_Rpassive(A1tpr, :), 20);
[mA1ra seA1ra] = jackmeanerr(A1_Ractive(A1tpr, :), 20);
[mA1tp seA1tp] = jackmeanerr(A1_Tpassive(A1tpr, :), 20);
[mA1ta seA1ta] = jackmeanerr(A1_Tactive(A1tpr, :), 20);

subplot(4,2,1);
shadedErrorBar(A1_t, mA1rp, seA1rp, '--b', 0);
hold on;
shadedErrorBar(A1_t, mA1ra, seA1ra, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('A1: Reference pass vs act');

subplot(4,2,2);
shadedErrorBar(A1_t, mA1tp, seA1tp, '--r', 0);
hold on;
shadedErrorBar(A1_t, mA1ta, seA1ta, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('A1: Target pass vs act');
text(1.2, yrange(2)*.8, ['N=' num2str(length(unique(A1_cells_used)))]);

% dPEG
if svd_subset == 0
  [mdPEGrp sedPEGrp] = jackmeanerr(dPEG_Rpassive(dPEGtpr, :), 20);
  [mdPEGra sedPEGra] = jackmeanerr(dPEG_Ractive (dPEGtpr, :), 20);
  [mdPEGtp sedPEGtp] = jackmeanerr(dPEG_Tpassive(dPEGtpr, :), 20);
  [mdPEGta sedPEGta] = jackmeanerr(dPEG_Tactive (dPEGtpr, :), 20);
  
  subplot(4,2,3);
  shadedErrorBar(dPEG_t, mdPEGrp, sedPEGrp, '--b', 0);
  hold on;
  shadedErrorBar(dPEG_t, mdPEGra, sedPEGra, 'b',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('dPEG: Reference pass vs act');
  
  subplot(4,2,4);
  shadedErrorBar(dPEG_t, mdPEGtp, sedPEGtp, '--r', 0);
  hold on;
  shadedErrorBar(dPEG_t, mdPEGta, sedPEGta, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('dPEG: Target pass vs act');
  text(1.2, yrange(2)*.8, ['N=' num2str(length(unique(dPEG_cells_used)))]);
end

% svd_dPEG
if svd_subset == 1
  
  [mdPEGrp sedPEGrp] = jackmeanerr(svd_dPEG_Rpassive(dPEGtpr, :), 20);
  [mdPEGra sedPEGra] = jackmeanerr(svd_dPEG_Ractive (dPEGtpr, :), 20);
  [mdPEGtp sedPEGtp] = jackmeanerr(svd_dPEG_Tpassive(dPEGtpr, :), 20);
  [mdPEGta sedPEGta] = jackmeanerr(svd_dPEG_Tactive (dPEGtpr, :), 20);  
  
  subplot(4,2,3);
  shadedErrorBar(svd_dPEG_t, mdPEGrp, sedPEGrp, '--b', 0);
  hold on;
  shadedErrorBar(svd_dPEG_t, mdPEGra, sedPEGra, 'b',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('dPEG: Reference pass vs act');
  
  subplot(4,2,4);
  shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, '--r', 0);
  hold on;
  shadedErrorBar(svd_dPEG_t, mdPEGta, sedPEGta, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('dPEG: Target pass vs act');
  text(1.2, yrange(2)*.8, ['N=' num2str(length(unique(dPEG_cells_used)))]);
end

% VPr
[mVPrrp seVPrrp] = jackmeanerr(VPr_Rpassive(VPrtpr, :), 20);
[mVPrra seVPrra] = jackmeanerr(VPr_Ractive (VPrtpr, :), 20);
[mVPrtp seVPrtp] = jackmeanerr(VPr_Tpassive(VPrtpr, :), 20);
[mVPrta seVPrta] = jackmeanerr(VPr_Tactive (VPrtpr, :), 20);

subplot(4,2,5);
shadedErrorBar(VPr_t, mVPrrp, seVPrrp, '--b', 0);
hold on;
shadedErrorBar(VPr_t, mVPrra, seVPrra, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 3.2]);
ylim(yrange);
title('VPr: Reference pass vs act');

subplot(4,2,6);
shadedErrorBar(VPr_t, mVPrtp, seVPrtp, '--r', 0);
hold on;
shadedErrorBar(VPr_t, mVPrta, seVPrta, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 3.2]);
ylim(yrange);
title('VPr: Target pass vs act');
text(2.2, yrange(2)*.8, ['N=' num2str(length(unique(VPr_cells_used)))]);

% dlFC
[mPFCrp sePFCrp] = jackmeanerr(PFC_Rpassive(PFCtpr, :), 20);
[mPFCra sePFCra] = jackmeanerr(PFC_Ractive (PFCtpr, :), 20);
[mPFCtp sePFCtp] = jackmeanerr(PFC_Tpassive(PFCtpr, :), 20);
[mPFCta sePFCta] = jackmeanerr(PFC_Tactive (PFCtpr, :), 20);

subplot(4,2,7);
shadedErrorBar(PFC_t, mPFCrp, sePFCrp, '--b', 0);
hold on;
shadedErrorBar(PFC_t, mPFCra, sePFCra, 'b',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('dlFC: Reference pass vs act');

subplot(4,2,8);
shadedErrorBar(PFC_t, mPFCtp, sePFCtp, '--r', 0);
hold on;
shadedErrorBar(PFC_t, mPFCta, sePFCta, 'r',   0);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('dlFC: Target pass vs act');
text(1.2, yrange(2)*.8, ['N=' num2str(length(unique(PFC_cells_used)))]);

suptitle('Target detection, passive vs. active');
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
hp1=shadedErrorBar(A1_t, mA1rp, seA1rp, '--b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1=shadedErrorBar(A1_t, mA1tp, seA1tp, '--r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('\bfPassive', 'FontSize', 14);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');
% set(gca, 'XTickLabel', []);

subplot(4,2,2);
hp1 = shadedErrorBar(A1_t, mA1ra, seA1ra, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(A1_t, mA1ta, seA1ta, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('\bfBehavior', 'FontSize', 14);
text(1.2, yrange(2)*.8, ['\bfN=' num2str(length(A1tpr))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);
% set(gca, 'XTickLabel', []);

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
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  % set(gca, 'XTickLabel', []);
  
  
  subplot(4,2,4);
  hp1 = shadedErrorBar(dPEG_t, mdPEGra, sedPEGra, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(dPEG_t, mdPEGta, sedPEGta, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  text(1.2, yrange(2)*.8, ['\bfN=' num2str(length(dPEGtpr))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);
  % set(gca, 'XTickLabel', []);
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
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  set(gca, 'box', 'off');
  ylabel('Normalized Firing Rate');
  % set(gca, 'XTickLabel', []);
  
  subplot(4,2,4);
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGra, sedPEGra, 'b', 0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  hold on;
  hp1 = shadedErrorBar(svd_dPEG_t, mdPEGta, sedPEGta, 'r',   0);
  set(hp1.edge, 'LineStyle', 'none');
  set(hp1.mainLine, 'LineWidth', 2);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
  line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  text(1.2, yrange(2)*.8, ['\bfN=' num2str(length(dPEGtpr))]);
  set(gca, 'box', 'off');
  set(gca, 'YTickLabel', []);
  % set(gca, 'XTickLabel', []);
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
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 3.2]);
ylim(yrange);
set(gca, 'box', 'off');
ylabel('Normalized Firing Rate');
% set(gca, 'XTickLabel', []);

subplot(4,2,6);
hp1 = shadedErrorBar(VPr_t, mVPrra, seVPrra, 'b', 0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
hold on;
hp1 = shadedErrorBar(VPr_t, mVPrta, seVPrta, 'r',   0);
set(hp1.edge, 'LineStyle', 'none');
set(hp1.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 3.2]);
ylim(yrange);
text(2.2, yrange(2)*.8, ['\bfN=' num2str(length(VPrtpr))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);
% set(gca, 'XTickLabel', []);

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
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1.8]);
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
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
text(1.2, yrange(2)*.8, ['\bfN=' num2str(length(PFCtpr))]);
set(gca, 'box', 'off');
set(gca, 'YTickLabel', []);
xlabel('Time (s)');

% suptitle('Target detection, reference vs. target');

set(gcf, 'PaperPositionMode', 'auto');
print('02', '-dpdf');


if latstat
  alfa = 0.05;
  twin = 36; % do test from bin 1 to 36 (1 s from onset)
  t0 = 7; % t=0 bin
  
  A1p_tar   = zeros(twin, 1);
  dPEGp_tar = zeros(twin, 1);
  VPrp_tar  = zeros(twin, 1);
  PFCp_tar  = zeros(twin, 1);
  
  A1p_ref   = zeros(twin, 1);
  dPEGp_ref = zeros(twin, 1);
  VPrp_ref  = zeros(twin, 1);
  PFCp_ref  = zeros(twin, 1);
  
  % target active vs target passive
  for n = 1 : twin
    A1p_tar(n)   = signrank(A1_Tpassive(:,n), A1_Tactive(:,n), 'alpha', alfa);
    dPEGp_tar(n) = signrank(dPEG_Tpassive(:,n), dPEG_Tactive(:,n), 'alpha', alfa);
    VPrp_tar(n)  = signrank(VPr_Tpassive(:,n), VPr_Tactive(:,n), 'alpha', alfa);
    PFCp_tar(n)  = signrank(PFC_Tpassive(:,n), PFC_Tactive(:,n), 'alpha', alfa);
  end
  
    % reference active vs reference passive
  for n = 1 : twin
    A1p_ref(n)   = signrank(A1_Rpassive(:,n), A1_Ractive(:,n), 'alpha', alfa);
    dPEGp_ref(n) = signrank(dPEG_Rpassive(:,n), dPEG_Ractive(:,n), 'alpha', alfa);
    VPrp_ref(n)  = signrank(VPr_Rpassive(:,n), VPr_Ractive(:,n), 'alpha', alfa);
    PFCp_ref(n)  = signrank(PFC_Rpassive(:,n), PFC_Ractive(:,n), 'alpha', alfa);
  end
  
  A1sig_ref   = findstr(A1p_ref(t0:end)', [1,1,1]);
  dPEGsig_ref = findstr(dPEGp_ref(t0:end)', [1,1,1]);
  VPrsig_ref  = findstr(VPrp_ref(t0:end)', [1,1,1]);
  PFCsig_ref  = findstr(PFCp_ref(t0:end)', [1,1,1]);
  
  A1sig_tar   = findstr(A1p_tar(t0:end)', [1,1,1]);
  dPEGsig_tar = findstr(dPEGp_tar(t0:end)', [1,1,1]);
  VPrsig_tar  = findstr(VPrp_tar(t0:end)', [1,1,1]);
  PFCsig_tar  = findstr(PFCp_tar(t0:end)', [1,1,1]);

keyboard;
[A1p_ref<(alfa/twin), dPEGp_ref<(alfa/twin), VPrp_ref<(alfa/twin), PFCp_ref< (alfa/twin)]
[A1p_tar<(alfa/twin), dPEGp_tar<(alfa/twin), VPrp_tar<(alfa/twin), PFCp_tar< (alfa/twin)]


[A1p_ref<alfa, dPEGp_ref<alfa, VPrp_ref<alfa, PFCp_ref< alfa]
[A1p_tar<alfa, dPEGp_tar<alfa, VPrp_tar<alfa, PFCp_tar< alfa]


[A1p_ref<alfa | A1p_tar<alfa, dPEGp_ref<alfa | dPEGp_tar<alfa, VPrp_ref<alfa | VPrp_tar<alfa, PFCp_ref< alfa | PFCp_tar< alfa]



end

if plot_ctx_figures
  %% plot 3: TORC vs FTC comparison with Passive
  
  
  figure;
  set(gcf, 'Position', [400 52 509 739]);
  
  
  % A1
  [mA1rr seA1rr] = jackmeanerr(A1_TOR(A1tpr, :), 20);
  [mA1tt seA1tt] = jackmeanerr(A1_FTC(A1tpr, :), 20);
  
  subplot(4,2,1);
  shadedErrorBar(A1_t, mA1rr, seA1rr, 'b', 0);
  hold on;
  shadedErrorBar(A1_t, mA1tt, seA1tt, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('A1: FTC vs TORC');
  
  subplot(4,2,2);
  shadedErrorBar(A1_t, mA1rp, seA1rp, 'b', 0);
  hold on;
  shadedErrorBar(A1_t, mA1tp, seA1tp, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('A1: Passive ref vs tar');
  text(1.2, yrange(2)*.8, ['N=' num2str(length(A1tpr))]);
  
  
  % dPEG
  if svd_subset == 0
    [mdPEGrr sedPEGrr] = jackmeanerr(dPEG_TOR(dPEGtpr, :), 20);
    [mdPEGtt sedPEGtt] = jackmeanerr(dPEG_FTC(dPEGtpr, :), 20);
    
    subplot(4,2,3);
    shadedErrorBar(dPEG_t, mdPEGrr, sedPEGrr, 'b', 0);
    hold on;
    shadedErrorBar(dPEG_t, mdPEGtt, sedPEGtt, 'r',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: FTC vs TORC');
    
    subplot(4,2,4);
    shadedErrorBar(dPEG_t, mdPEGrp, sedPEGrp, 'b', 0);
    hold on;
    shadedErrorBar(dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: Passive ref vs tar');
    text(1.2, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % svd_dPEG
  if svd_subset == 1
    [mdPEGrr sedPEGrr] = jackmeanerr(svd_dPEG_TOR(dPEGtpr, :), 20);
    [mdPEGtt sedPEGtt] = jackmeanerr(svd_dPEG_FTC(dPEGtpr, :), 20);
    
    subplot(4,2,3);
    shadedErrorBar(svd_dPEG_t, mdPEGrr, sedPEGrr, 'b', 0);
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtt, sedPEGtt, 'r',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: FTC vs TORC');
    
    subplot(4,2,4);
    shadedErrorBar(svd_dPEG_t, mdPEGrp, sedPEGrp, 'b', 0);
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: Passive ref vs tar');
    text(1.2, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % VPr
  [mVPrrr seVPrrr] = jackmeanerr(VPr_TOR(VPrtpr, :), 20);
  [mVPrtt seVPrtt] = jackmeanerr(VPr_FTC(VPrtpr, :), 20);
  
  subplot(4,2,5);
  shadedErrorBar(VPr_t, mVPrrr, seVPrrr, 'b', 0);
  hold on;
  shadedErrorBar(VPr_t, mVPrtt, seVPrtt, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('VPr: FTC vs TORC');
  
  subplot(4,2,6);
  shadedErrorBar(VPr_t, mVPrrp, seVPrrp, 'b', 0);
  hold on;
  shadedErrorBar(VPr_t, mVPrtp, seVPrtp, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('VPr: Passive ref vs tar');
  text(1.2, yrange(2)*.8, ['N=' num2str(length(VPrtpr))]);
  
  
  % PFC
  [mPFCrr sePFCrr] = jackmeanerr(PFC_TOR(PFCtpr, :), 20);
  [mPFCtt sePFCtt] = jackmeanerr(PFC_FTC(PFCtpr, :), 20);
  
  subplot(4,2,7);
  shadedErrorBar(PFC_t, mPFCrr, sePFCrr, 'b', 0);
  hold on;
  shadedErrorBar(PFC_t, mPFCtt, sePFCtt, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('PFC: FTC vs TORC');
  
  subplot(4,2,8);
  shadedErrorBar(PFC_t, mPFCrp, sePFCrp, 'b', 0);
  hold on;
  shadedErrorBar(PFC_t, mPFCtp, sePFCtp, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  xlim([-0.2 1.8]);
  ylim(yrange);
  hold off;
  text(1.2, yrange(2)*.8, ['N=' num2str(length(PFCtpr))]);
  
  % suptitle('TORC vs. tone, out and in task context');
  
  set(gcf, 'PaperPositionMode', 'auto');
  print('03', '-dpdf');
  
  %% plot 4: TORC vs pass ref / FTC vs pass tar
  figure;
  set(gcf, 'Position', [500 52 509 739]);
  
  % A1
  subplot(4,2,1);
  shadedErrorBar(A1_t, mA1rr, seA1rr, '--b', 0);
  hold on;
  shadedErrorBar(A1_t, mA1rp, seA1rp, 'b',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('A1: TORC vs pass TORC');
  
  subplot(4,2,2);
  shadedErrorBar(A1_t, mA1tt, seA1tt, '--r', 0);
  hold on;
  shadedErrorBar(A1_t, mA1tp, seA1tp, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('A1: FTC vs pass target');
  text(1.2, yrange(2)*.8, ['N=' num2str(length(A1tpr))]);
  
  
  % dPEG
  if svd_subset == 0
    subplot(4,2,3);
    shadedErrorBar(dPEG_t, mdPEGrr, sedPEGrr, '--b', 0);
    hold on;
    shadedErrorBar(dPEG_t, mdPEGrp, sedPEGrp, 'b',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: TORC vs pass TORC');
    
    subplot(4,2,4);
    shadedErrorBar(dPEG_t, mdPEGtt, sedPEGtt, '--r', 0);
    hold on;
    shadedErrorBar(dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: FTC vs pass target');
    text(1.2, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % svd_dPEG
  if svd_subset == 1
    subplot(4,2,3);
    shadedErrorBar(svd_dPEG_t, mdPEGrr, sedPEGrr, '--b', 0);
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, 'b',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: TORC vs pass TORC');
    
    subplot(4,2,4);
    shadedErrorBar(svd_dPEG_t, mdPEGtt, sedPEGtt, '--r', 0);
    hold on;
    shadedErrorBar(svd_dPEG_t, mdPEGtp, sedPEGtp, 'r',   0);
    line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
    hold off;
    xlim([-0.2 1.8]);
    ylim(yrange);
    title('dPEG: FTC vs pass target');
    text(1.2, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);
  end
  
  % VPr
  subplot(4,2,5);
  shadedErrorBar(VPr_t, mVPrrr, seVPrrr, '--b', 0);
  hold on;
  shadedErrorBar(VPr_t, mVPrrp, seVPrrp, 'b',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('VPr: TORC vs pass TORC');
  
  subplot(4,2,6);
  shadedErrorBar(VPr_t, mVPrtt, seVPrtt, '--r', 0);
  hold on;
  shadedErrorBar(VPr_t, mVPrtp, seVPrtp, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('VPr: FTC vs pass target');
  text(1.2, yrange(2)*.8, ['N=' num2str(length(VPrtpr))]);
  
  
  % PFC
  subplot(4,2,7);
  shadedErrorBar(PFC_t, mPFCrr, sePFCrr, '--b', 0);
  hold on;
  shadedErrorBar(PFC_t, mPFCrp, sePFCrp, 'b',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  hold off;
  xlim([-0.2 1.8]);
  ylim(yrange);
  title('PFC: TORC vs pass TORC');
  
  subplot(4,2,8);
  shadedErrorBar(PFC_t, mPFCtt, sePFCtt, '--r', 0);
  hold on;
  shadedErrorBar(PFC_t, mPFCtp, sePFCtp, 'r',   0);
  line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
  xlim([-0.2 1.8]);
  ylim(yrange);
  hold off;
  title('PFC: FTC vs pass target');
  text(1.2, yrange(2)*.8, ['N=' num2str(length(PFCtpr))]);
  
  suptitle('Out of task vs. passive task');
  
  set(gcf, 'PaperPositionMode', 'auto');
  print('04', '-dpdf');
end

%% Plot 5: Naive vs Trained figure
if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.8 1.3];
elseif ~popnorm && ~normalize_each
  yrange = [-20 20];
end

figure;
set(gcf, 'Position', [21 48 693 741]);
set(gcf, 'Color', 'w');
% A1
[mnA1rp senA1rp] = jackmeanerr(nA1_Rpassive, 20);
[mnA1tp senA1tp] = jackmeanerr(nA1_Tpassive, 20);

subplot(3,3,1);
h5 = shadedErrorBar(nA1_t, mnA1rp, senA1rp, 'b', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
hold on;
h5 = shadedErrorBar(nA1_t, mnA1tp, senA1tp, 'r', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
xlim([-0.2 1]);
ylim(yrange);
hold off;
title('\bfNaive Tone Detection');
text(0.8, yrange(2)*.8, ['\bfN=' num2str(size(nA1_Tpassive,1))]);
set(gca, 'box', 'off');
ylabel('Normalized firing rate');
text(-0.8, 0.4, '\bfA1', 'FontSize', 12);
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
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1]);
ylim(yrange);
title('\bfTrained Passive Tone Detection');
set(gca, 'box', 'off');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Layer', 'top');

subplot(3,3,3);
h5 = shadedErrorBar(A1_t, mA1ra, seA1ra, 'b', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
hold on;
h5 = shadedErrorBar(A1_t, mA1ta, seA1ta, 'r',   0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1]);
ylim(yrange);
ylim(yrange);
title('\bfBehaving Tone Detection');
text(0.7, yrange(2)*.8, ['\bfN=' num2str(length(A1tpr))]);
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

[mndPEGrp sendPEGrp] = jackmeanerr(ndPEG_Rpassive, 20);
[mndPEGtp sendPEGtp] = jackmeanerr(ndPEG_Tpassive, 20);

subplot(3,3,4);
h5 = shadedErrorBar(ndPEG_t, mndPEGrp, sendPEGrp, 'b', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
hold on;
h5 = shadedErrorBar(ndPEG_t, mndPEGtp, sendPEGtp, 'r', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
xlim([-0.2 1]);
ylim(yrange);
hold off;
text(0.8, yrange(2)*.8, ['\bfN=' num2str(size(ndPEG_Tpassive,1))]);
set(gca, 'box', 'off');
ylabel('Normalized firing rate');
text(-0.9, 0.4, '\bfdPEG', 'FontSize', 12);
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
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1]);
ylim(yrange);
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
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1]);
ylim(yrange);
set(gca, 'box', 'off');
text(0.7, yrange(2)*.8, ['\bfN=' num2str(length(dPEGtpr))]);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Layer', 'top');

% VPr
[mnVPrrp senVPrrp] = jackmeanerr(nVPr_Rpassive, 20);
[mnVPrtp senVPrtp] = jackmeanerr(nVPr_Tpassive, 20);

subplot(3,3,7);
h5 = shadedErrorBar(nVPr_t(1:60), mnVPrrp(1:60), senVPrrp(1:60), 'b', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
hold on;
h5 = shadedErrorBar(nVPr_t(1:60), mnVPrtp(1:60), senVPrtp(1:60), 'r', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
xlim([-0.2 1]);
ylim(yrange);
hold off;
text(0.8, yrange(2)*.8, ['\bfN=' num2str(size(nVPr_Tpassive,1))]);
set(gca, 'box', 'off');
set(gca, 'Layer', 'top');
xlabel('Time (s)');
ylabel('Normalized firing rate');
text(-0.9, 0.4, '\bfVPr', 'FontSize', 12);

subplot(3,3,8);
VPcorr = max([mVPrrp(1:60) mVPrra(1:60) mVPrtp(1:60) mVPrta(1:60)]);
h5 = shadedErrorBar(VPr_t(1:60), mVPrrp(1:60)./VPcorr, seVPrrp(1:60), 'b',   0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
hold on;
h5 = shadedErrorBar(VPr_t(1:60), mVPrtp(1:60)./VPcorr, seVPrtp(1:60), 'r',   0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1]);
ylim(yrange);
set(gca, 'box', 'off');
xlabel('Time (s)');
set(gca, 'YTickLabel', []);
set(gca, 'Layer', 'top');

subplot(3,3,9);
h5 = shadedErrorBar(VPr_t(1:60), mVPrra(1:60)./VPcorr, seVPrra(1:60), 'b', 0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
hold on;
h5 = shadedErrorBar(VPr_t(1:60), mVPrta(1:60)./VPcorr, seVPrta(1:60), 'r',   0);
set(h5.edge, 'LineStyle', 'none');
set(h5.mainLine, 'LineWidth', 2);
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1]);
ylim(yrange);
set(gca, 'box', 'off');
text(0.7, yrange(2)*.8, ['\bfN=' num2str(length(VPrtpr))]);
xlabel('Time (s)');
set(gca, 'YTickLabel', []);
set(gca, 'Layer', 'top');

set(gcf, 'PaperPositionMode', 'auto');
print('05', '-dpdf');

%% trying to do stats on naive vs trained-passive
% nDtr = [nA1_Dtrp; ndPEG_Dtrp; nVPr_Dtrp];
% Dtr = [A1_Dtrp; dPEG_Dtrp; VPr_Dtrp];
% Dtrg= [ones(length(A1_Dtrp),1); ones(length(dPEG_Dtrp),1).*2; ones(length(VPr_Dtrp),1).*3];
% nDtrg= [ones(length(nA1_Dtrp),1); ones(length(ndPEG_Dtrp),1).*2; ones(length(nVPr_Dtrp),1).*3];
% Dtr_area = [nDtrg;Dtrg];
% Dtr_all = [nDtr;Dtr];
% Dtr_tr   = [ones(length(nDtr),1); ones(length(Dtr),1).*2];
% 
% [p tbl stats] = anovan(Dtr_all,{Dtr_area, Dtr_tr}, 'model' ,2, 'varnames', {'Area','Training'});
% figure;c=multcompare(stats, 'Dimension', 1)
% figure;multcompare(stats, 'Dimension', 2)
% figure;c=multcompare(stats, 'Dimension', [1 2])
% TD
%   Source          Sum Sq.   d.f.   Mean Sq.     F     Prob>F
% ------------------------------------------------------------
%   Area             26.615     2    13.3077    11.54   0     
%   Training          7.821     1     7.8214     6.78   0.0094
%   Area*Training     0.335     2     0.1675     0.15   0.8648
%   Error           846.672   734     1.1535                  
%   Total           894.834   739                             

% CD
%   Source          Sum Sq.   d.f.   Mean Sq.    F     Prob>F
% -----------------------------------------------------------
%   Area              0.259     2    0.12949    0.11   0.8983
%   Training          7.035     1    7.03537    5.83   0.016 
%   Area*Training     3.254     2    1.62709    1.35   0.2603
%   Error           885.625   734    1.20657                 
%   Total           902.724   739                            


%% Plot 6: Tar-Ref differences

if only_pos_tar == -1
  yrange([-0.5 0.5]);
else
  yrange =[-0.2 0.6];
end

if popnorm
  yrange = [-0.8 1.3];
elseif ~popnorm && ~normalize_each
  yrange = [-20 20];
end

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

ntempdiff = nA1_Tpassive - nA1_Rpassive; % naive
[nA1_Dtr_passive  nA1_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
nA1_Dtrp = nanmean(ntempdiff(:,scatwin),2);

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
  
  ntempdiff = ndPEG_Tpassive - ndPEG_Rpassive; % naive
  [ndPEG_Dtr_passive  ndPEG_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
  ndPEG_Dtrp = nanmean(ntempdiff(:,scatwin),2);
  
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
  
  ntempdiff = ndPEG_Tpassive - ndPEG_Rpassive; % naive
  [ndPEG_Dtr_passive  ndPEG_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
  ndPEG_Dtrp = nanmean(ntempdiff(:,scatwin),2);
  
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

ntempdiff = nVPr_Tpassive - nVPr_Rpassive; % naive
[nVPr_Dtr_passive  nVPr_Dtr_passive_sem] = jackmeanerr(ntempdiff, 20);
nVPr_Dtrp = nanmean(ntempdiff(:,scatwin),2);

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
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
hold off;
title('A1 Change in reference response');

subplot(4,3,2);
shadedErrorBar(A1_t, A1_Dap_tar, A1_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Change in target response');

subplot(4,3,3);
shadedErrorBar(A1_t, mA1_tarenhancement, seA1_tarenhancement, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Tar change - Ref change');
text(1.2, yrange(2)*.8, ['N=' num2str(length(A1tpr))]);


% dPEG
subplot(4,3,4);
shadedErrorBar(dPEG_t, dPEG_Dap_ref, dPEG_Dap_ref_sem, 'b', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('dPEG Change in reference response');

subplot(4,3,5);
shadedErrorBar(dPEG_t, dPEG_Dap_tar, dPEG_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Change in target response');

subplot(4,3,6);
shadedErrorBar(dPEG_t, mdPEG_tarenhancement, sedPEG_tarenhancement, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Tar change - Ref change');
text(1.2, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);


% VPr
subplot(4,3,7);
shadedErrorBar(VPr_t, VPr_Dap_ref, VPr_Dap_ref_sem, 'b', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
ylim(yrange);
xlim([-0.2 3.2]);
title('VPr Change in reference response');

subplot(4,3,8);
shadedErrorBar(VPr_t, VPr_Dap_tar, VPr_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
ylim(yrange);
xlim([-0.2 3.2]);
title('Change in target response');

subplot(4,3,9);
shadedErrorBar(VPr_t, mVPr_tarenhancement, seVPr_tarenhancement, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2 2],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([2.4 2.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([2.8 2.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 3.2]);
ylim(yrange);
title('Tar change - Ref change');
text(2.2, yrange(2)*.8, ['N=' num2str(length(VPrtpr))]);


% PFC
subplot(4,3,10);
shadedErrorBar(PFC_t, PFC_Dap_ref, PFC_Dap_ref_sem, 'b', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
title('PFC Change in reference response');

subplot(4,3,11);
shadedErrorBar(PFC_t, PFC_Dap_tar, PFC_Dap_tar_sem, 'r',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
title('Change in target response');

subplot(4,3,12);
shadedErrorBar(PFC_t, mPFC_tarenhancement, sePFC_tarenhancement, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Tar change - Ref change');
text(1.2, yrange(2)*.8, ['N=' num2str(length(PFCtpr))]);

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
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
hold off;
title('A1 Passive tar-ref contrast');

subplot(4,3,2);
shadedErrorBar(A1_t, A1_Dtr_active, A1_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
title('Active tar-ref contrast');

subplot(4,3,3);
shadedErrorBar(A1_t, mA1_contrastenh, seA1_contrastenh, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Tar - Ref contrast change');
text(1.2, yrange(2)*.8, ['N=' num2str(length(A1tpr))]);


% dPEG
subplot(4,3,4);
shadedErrorBar(dPEG_t, dPEG_Dtr_passive, dPEG_Dtr_passive_sem, '--m', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
hold off;
title('dPEG Passive tar-ref contrast');

subplot(4,3,5);
shadedErrorBar(dPEG_t, dPEG_Dtr_active, dPEG_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
title('Active tar-ref contrast');

subplot(4,3,6);
shadedErrorBar(dPEG_t, mdPEG_contrastenh, sedPEG_contrastenh, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Tar - Ref contrast change');
text(1.2, yrange(2)*.8, ['N=' num2str(length(dPEGtpr))]);


% VPr
subplot(4,3,7);
shadedErrorBar(VPr_t, VPr_Dtr_passive, VPr_Dtr_passive_sem, '--m', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
hold off;
title('VPr Passive tar-ref contrast');

subplot(4,3,8);
shadedErrorBar(VPr_t, VPr_Dtr_active, VPr_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
title('Active tar-ref contrast');

subplot(4,3,9);
shadedErrorBar(VPr_t, mVPr_contrastenh, seVPr_contrastenh, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Tar - Ref contrast change');
text(1.2, yrange(2)*.8, ['N=' num2str(length(VPrtpr))]);


% PFC
subplot(4,3,10);
shadedErrorBar(PFC_t, PFC_Dtr_passive, PFC_Dtr_passive_sem, '--m', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
title('PFC Passive tar-ref contrast');

subplot(4,3,11);
shadedErrorBar(PFC_t, PFC_Dtr_active, PFC_Dtr_active_sem, 'm',   0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
ylim(yrange);
xlim([-0.2 1.8]);
title('Active tar-ref contrast');

subplot(4,3,12);
shadedErrorBar(PFC_t, mPFC_contrastenh, sePFC_contrastenh, 'k', 0);
hold on;
line([0 0],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1 1],yrange, 'LineStyle', '--', 'Color', [.5 .5 .5]);
line([1.4 1.4],yrange, 'LineStyle', '--', 'Color', 'g');
line([1.8 1.8],yrange, 'LineStyle', '--', 'Color', 'g');
hold off;
xlim([-0.2 1.8]);
ylim(yrange);
title('Tar - Ref contrast change');
text(1.2, yrange(2)*.8, ['N=' num2str(length(PFCtpr))]);

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
elseif ~popnorm && ~normalize_each
  yrange = [-20 20];
end

figure;
subplot(2,2,1);
scatter(A1_Dtrp, A1_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange)
ylim(yrange)
xlabel('Passive target-reference')
ylabel('Active target-reference')
title('A1')
hold on
line(yrange, yrange, 'Color', 'k')
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k')
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k')

subplot(2,2,2);
scatter(dPEG_Dtrp, dPEG_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange)
ylim(yrange)
xlabel('Passive target-reference')
ylabel('Active target-reference')
title('dPEG')
hold on
line(yrange, yrange, 'Color', 'k')
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k')
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k')

subplot(2,2,3);
scatter(VPr_Dtrp, VPr_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange)
ylim(yrange)
xlabel('Passive target-reference')
ylabel('Active target-reference')
title('VPr')
hold on
line(yrange, yrange, 'Color', 'k')
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k')
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k')

subplot(2,2,4);
scatter(PFC_Dtrp, PFC_Dtra, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
xlim(yrange)
ylim(yrange)
xlabel('Passive target-reference')
ylabel('Active target-reference')
title('PFC')
hold on
line(yrange, yrange, 'Color', 'k')
line([0 0], yrange, 'LineStyle', '--', 'Color', 'k')
line(yrange, [0 0], 'LineStyle', '--', 'Color', 'k')

set(gcf, 'PaperPositionMode', 'auto');
print('08', '-dpdf');

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
print('09', '-dpdf');

%% figure 10 : naive, passive and active tar-ref contrast histograms

lefthlf = [1: floor(length(edges)/2)];
righthlf = [floor(length(edges)/2)+1 : length(edges)];
leftclr = [64/255, 186/255, 64/255];
rightclr = [0/255,127/255,0/255];

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
print('10', '-dpdf');


%% Figure 11 definitions
%  Analysis of post-target response
% ptwin1= [39:39+10];%60]; %1.08 to 1.4167s
% ptwin2= [69:69+10];%90]; %2.0833 to 2.4167s
% ptwin1= [39:57]; %1.08 to 1.7833s
% ptwin2= [69:87]; %2.0833 to 2.7833s
% ptwin1= [37:58]; %1.0167 to 1.7167s
% ptwin2= [67:88]; %2.0167 to 2.7167s
ptwin1= [38:60]; %1.05 to 1.78s
ptwin2= [68:90]; %2.05 to 2.78s
if popnorm == 1
  xrange = [-3 3];
elseif popnorm == 0 && normalize_each == 0
  xrange = [-20 20];
elseif popnorm == 0 && normalize_each == 1
  xrange = [-1 1];
end


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

%edges = [-1:0.1:1];
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
print('11', '-dpdf');

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
multcompare(pt_stats);


%% figure 12 : change in reference vs change in target scatter plot
if normalize_each == 0
  xrange = [-10 10];
end
if popnorm == 1
  xrange = [-5 5];
else
  xrange = [-15 15];
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
% suptitle('Spontaneous activity - Tone detection');
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print('13', '-dpdf');

return;

%% organize and save variables for PTD - CLK comparisons

PTD_A1_cells_used = A1_cells_used;
PTD_A1_TEI = A1_Dtra - A1_Dtrp;
PTD_A1_Dtrp = A1_Dtrp;
PTD_A1_Dtra = A1_Dtra;
PTD_A1_at = A1_at;
PTD_A1_Dapr = A1_Dapr;
PTD_A1_Dapt = A1_Dapt;

PTD_dPEG_cells_used = dPEG_cells_used;
PTD_dPEG_TEI = dPEG_Dtra - dPEG_Dtrp;
PTD_dPEG_Dtrp = dPEG_Dtrp;
PTD_dPEG_Dtra = dPEG_Dtra;
PTD_dPEG_at = dPEG_at;
PTD_dPEG_Dapr = dPEG_Dapr;
PTD_dPEG_Dapt = dPEG_Dapt;

PTD_VPr_cells_used = VPr_cells_used;
PTD_VPr_TEI = VPr_Dtra - VPr_Dtrp;
PTD_VPr_Dtrp = VPr_Dtrp;
PTD_VPr_Dtra = VPr_Dtra;
PTD_VPr_at = VPr_at;
PTD_VPr_Dapr = VPr_Dapr;
PTD_VPr_Dapt = VPr_Dapt;

PTD_PFC_cells_used = PFC_cells_used;
PTD_PFC_TEI = PFC_Dtra - PFC_Dtrp;
PTD_PFC_Dtrp = PFC_Dtrp;
PTD_PFC_Dtra = PFC_Dtra;
PTD_PFC_at = PFC_at;
PTD_PFC_Dapr = PFC_Dapr;
PTD_PFC_Dapt = PFC_Dapt;

save('PTD_scatvals.mat', 'PTD_*');



%% figure XX : histograms comparing tar-ref contrast in naive, passive and active
% 
% 
% figure;
% set(gcf, 'Position', [2 77 559 747]);
% subplot(4,1,1);
% nA1_h = histc(nA1_Dtrp, edges);
% pA1_h = histc(A1_Dtrp, edges);
% aA1_h = histc(A1_Dtra, edges);
% maxyval = max(max([nA1_h pA1_h aA1_h]))*1.2;
% hold on;
% plot(edges, nA1_h, 'k', 'LineWidth', 3);
% plot(edges, pA1_h, 'b', 'LineWidth', 3);
% plot(edges, aA1_h, 'r', 'LineWidth', 3);
% line([0 0], [0 maxyval], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
% hold off;
% legend('Naive', 'Passive', 'Active');
% xlim(xrange);
% ylim([0 maxyval]);
% ylabel('Number of cells');
% title('A1');
% 
% subplot(4,1,2);
% ndPEG_h = histc(ndPEG_Dtrp, edges);
% pdPEG_h = histc(dPEG_Dtrp, edges);
% adPEG_h = histc(dPEG_Dtra, edges);
% maxyval = max(max([ndPEG_h pdPEG_h adPEG_h]))*1.2;
% hold on;
% plot(edges, ndPEG_h, 'k', 'LineWidth', 3);
% plot(edges, pdPEG_h, 'b', 'LineWidth', 3);
% plot(edges, adPEG_h, 'r', 'LineWidth', 3);
% line([0 0], [0 maxyval], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
% hold off;
% legend('Naive', 'Passive', 'Active');
% xlim(xrange);
% ylim([0 maxyval]);
% ylabel('Number of cells');
% title('dPEG');
% 
% subplot(4,1,3);
% nVPr_h = histc(nVPr_Dtrp, edges);
% pVPr_h = histc(VPr_Dtrp, edges);
% aVPr_h = histc(VPr_Dtra, edges);
% maxyval = max(max([nVPr_h pVPr_h aVPr_h]))*1.2;
% hold on;
% plot(edges, nVPr_h, 'k', 'LineWidth', 3);
% plot(edges, pVPr_h, 'b', 'LineWidth', 3);
% plot(edges, aVPr_h, 'r', 'LineWidth', 3);
% line([0 0], [0 maxyval], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
% hold off;
% legend('Naive', 'Passive', 'Active');
% xlim(xrange);
% ylim([0 maxyval]);
% ylabel('Number of cells');
% title('VPr');
% 
% subplot(4,1,4);
% %nPFC_h = histc(nPFC_Dtrp, edges);
% pPFC_h = histc(PFC_Dtrp, edges);
% aPFC_h = histc(PFC_Dtra, edges);
% maxyval = max(max([pPFC_h aPFC_h]))*1.2;
% hold on;
% %plot(edges, nPFC_h, 'k', 'LineWidth', 3);
% plot(edges, pPFC_h, 'b', 'LineWidth', 3);
% plot(edges, aPFC_h, 'r', 'LineWidth', 3);
% line([0 0], [0 maxyval], 'Color', 'k', 'LineWidth' , 2, 'LineStyle', '--');
% hold off;
% legend('Passive', 'Active');
% xlim(xrange);
% ylim([0 maxyval]);
% ylabel('Number of cells');
% title('PFC');


%% Notas reunión:
% conteo de neuronas totales,
% con respuestas pasivas
% con efectos
% sin efectos
% sin respuestas pasivas
% con efectos
% sin efectos
% PERO SOLO ENFATIZAR dlFC y VPr, porque en dPEG y A1 no necesariamente se
% registro de neuronas sin tuning/respuestas pasivas (sesgo del resultado)

% contar (y plotear) las neuronas con supresión (actualmente uso solo
% respuestas positivas a target).

% figura "target detection, reference vs. target" agregarle: todas las
% neuronas (positivas, y negativas) y otra columna con todas las neuronas
% suprimidas.
% Además, agregarle columnas de delta passive vs. delta active y de delta
% reference vs delta target. Tambien una curva de diferencia de
% diferencias.

% componer figuras de naive vs trained: 3 cols, naive passive, trained
% passive, trained active.

% Hacer tambien esto con CLKs

% try compare evoked and spont firing rate in all areas



%% try to correlate behavioral performance with off response

runclass = 'PTD';
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


%% plot results of post-target response vs behavior performance
behavior_measure = DR;

figure;scatter(post_target_index, behavior_measure)
%xlim([-5 5]);
ylim([0 100]);
line([0 0], [0 100], 'LineStyle', '--', 'Color', 'k');

figure; bar([median(post_target_index(behavior_measure<25)), median(post_target_index(behavior_measure>25 & behavior_measure<50)), median(post_target_index(behavior_measure>50 & behavior_measure<75)), median(post_target_index(behavior_measure>75))])

groups = zeros(length(post_target_index), 1);
groups(behavior_measure<25)                        = 1;
groups(behavior_measure>=25 & behavior_measure<50) = 2;
groups(behavior_measure>=50 & behavior_measure<75) = 3;
groups(behavior_measure>=75)                       = 4;



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


%%
figure; bar([median(post_target_index(behavior_measure<50)),  median(post_target_index(behavior_measure>=50))])

figure; 
chanced= histc(post_target_index(groups==1), [-6:.5:6])
plot([-6:0.5:6], chanced./max(chanced)); hold on


goodd= histc(post_target_index(groups==3), [-6:1:6])
plot([-6:1:6], goodd./max(goodd), 'r');


%% TRY to come up with an index of change
% Jonathan's idea of calculating the ratio if the change in the area under the curves
% of figure 4 in the ms

% during 0.1 to 0.45 from onset (scatwin)
anawin = [7:37];
A1enh = (trapz(A1_t(anawin), mA1ta(anawin)) - trapz(A1_t(anawin), mA1ra(anawin))) / (trapz(A1_t(anawin), mA1tp(anawin)) - trapz(A1_t(anawin), mA1rp(anawin)));
dPEGenh = (trapz(dPEG_t(anawin), mdPEGta(anawin)) - trapz(dPEG_t(anawin), mdPEGra(anawin))) / (trapz(dPEG_t(anawin), mdPEGtp(anawin)) - trapz(dPEG_t(anawin), mdPEGrp(anawin)));
VPrenh = (trapz(VPr_t(anawin), mVPrta(anawin)) - trapz(VPr_t(anawin), mVPrra(anawin))) / (trapz(VPr_t(anawin), mVPrtp(anawin)) - trapz(VPr_t(anawin), mVPrrp(anawin)));
PFCenh = (trapz(PFC_t(anawin), mPFCta(anawin)) - trapz(PFC_t(anawin), mPFCra(anawin))) / (trapz(PFC_t(anawin), mPFCtp(anawin)) - trapz(PFC_t(anawin), mPFCrp(anawin)));
figure;
plot([A1enh, dPEGenh, VPrenh, PFCenh])
hold on;
plot([A1enh, dPEGenh, VPrenh, PFCenh], 'ro');


% from onset to offset+.4s, separating passive and active
anawin  = [7:48];
anawin1 = [7:78];
lastt   = A1_t(anawin(end));
lastt1  = VPr_t(anawin1(end));

A1enhP = trapz(A1_t(anawin), mA1tp(anawin)) - trapz(A1_t(anawin), mA1rp(anawin));
A1enhA = trapz(A1_t(anawin), mA1ta(anawin)) - trapz(A1_t(anawin), mA1ra(anawin));
dPEGenhP = trapz(dPEG_t(anawin), mdPEGtp(anawin)) - trapz(dPEG_t(anawin), mdPEGrp(anawin));
dPEGenhA = trapz(dPEG_t(anawin), mdPEGta(anawin)) - trapz(dPEG_t(anawin), mdPEGra(anawin));
VPrenhP = trapz(VPr_t(anawin1), mVPrtp(anawin1)) - trapz(VPr_t(anawin1), mVPrrp(anawin1));
VPrenhA = trapz(VPr_t(anawin1), mVPrta(anawin1)) - trapz(VPr_t(anawin1), mVPrra(anawin1));
PFCenhP = trapz(PFC_t(anawin), mPFCtp(anawin)) - trapz(PFC_t(anawin), mPFCrp(anawin));
PFCenhA = trapz(PFC_t(anawin), mPFCta(anawin)) - trapz(PFC_t(anawin), mPFCra(anawin));

A1enhP   = A1enhP   / lastt;
A1enhA   = A1enhA   / lastt;
dPEGenhP = dPEGenhP / lastt;
dPEGenhA = dPEGenhA / lastt;
VPrenhP  = VPrenhP  / lastt1;
VPrenhA  = VPrenhA  / lastt1;
PFCenhP  = PFCenhP  / lastt;
PFCenhA  = PFCenhA  / lastt;

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

% from onset to offset+.4s, ref passive vs target active
anawin  = [7:48];
anawin1 = [7:78];
lastt   = A1_t(anawin(end));
lastt1  = VPr_t(anawin1(end));

A1enh = trapz(A1_t(anawin), mA1ta(anawin)) - trapz(A1_t(anawin), mA1rp(anawin));
dPEGenh = trapz(dPEG_t(anawin), mdPEGta(anawin)) - trapz(dPEG_t(anawin), mdPEGrp(anawin));
VPrenh = trapz(VPr_t(anawin1), mVPrta(anawin1)) - trapz(VPr_t(anawin1), mVPrrp(anawin1));
PFCenh = trapz(PFC_t(anawin), mPFCta(anawin)) - trapz(PFC_t(anawin), mPFCrp(anawin));

A1enh   = A1enh   / lastt;
dPEGenh = dPEGenh / lastt;
VPrenh  = VPrenh  / lastt1;
PFCenh  = PFCenh  / lastt;

figure;
bar([A1enh, dPEGenh, VPrenh, PFCenh]);

% % instead of helping with the previous attempt to find an index, PY and DD
% % suggested measuring the mean variance in TORC responses...
% 
% A1pv = nanvar(A1_Rpassive,0, 2);
% A1av = nanvar(A1_Ractive, 0, 2);
% dPEGpv = nanvar(dPEG_Rpassive,0, 2);
% dPEGav = nanvar(dPEG_Ractive, 0, 2);
% VPrpv = nanvar(VPr_Rpassive,0, 2);
% VPrav = nanvar(VPr_Ractive, 0, 2);
% PFCpv = nanvar(PFC_Rpassive,0, 2);
% PFCav = nanvar(PFC_Ractive, 0, 2);
% 
% [nanmean(A1pv) nanmean(A1av) nanmean(dPEGpv) nanmean(dPEGav) nanmean(VPrpv) nanmean(VPrav) nanmean(PFCpv) nanmean(PFCav)]
% 
% % but this should not work; we need the variance trial by trial, on the
% % hypothesis that for VP and FC the TORCs become being treated as one type
% % of sound, while A1 does code for the individual TORC properties.


%% # neurons with no response in passive, with effects

A1pfcl   = (sum(~A1_pass_mods)/length(A1_pass_mods)) * 100;
dPEGpfcl = (sum(~dPEG_pass_mods)/length(dPEG_pass_mods)) * 100;
VPrpfcl  = (sum(~VPr_pass_mods)/length(VPr_pass_mods)) * 100;
PFCpfcl  = (sum(~PFC_pass_mods)/length(PFC_pass_mods)) * 100;

figure('name', 'PFC-like neurons per Area');
bar([1,2,3,4], [A1pfcl, dPEGpfcl, VPrpfcl, PFCpfcl]);

%% separate population in first and last third to see if enhanced contrast appears with training
wascells = NaN(length(VPr_cells_used), 1);
safcells = NaN(length(VPr_cells_used), 1);
avocells = NaN(length(VPr_cells_used), 1);
lemcells = NaN(length(VPr_cells_used), 1);

for n = 1 : length(VPr_cells_used)
  if strcmp(VPr_cells_used{n}(1:3), 'was')
    wascells(n) = n;
  elseif strcmp(VPr_cells_used{n}(1:3), 'saf')
    safcells(n) = n;
  elseif strcmp(VPr_cells_used{n}(1:3), 'avo')
    avocells(n) = n;
  elseif strcmp(VPr_cells_used{n}(1:3), 'lem')
    lemcells(n) = n;
  end
end

wascells(isnan(wascells))=[];
safcells(isnan(safcells))=[];
lemcells(isnan(lemcells))=[];
avocells(isnan(avocells))=[];

figure; hold on;
listchunk = floor(length(safcells) / 10);
lvl = 0;
for i = 1: listchunk:length(safcells)
  stop = i+(listchunk-1);
  if stop > length (safcells),
    stop = length(safcells)
  elseif stop == length(safcells)
    break;
  end
  plot(VPr_t, nanmean(VPr_Rpassive(safcells(i:stop),:))+lvl); hold on;
  lvl = lvl+1;
end
title('saffon')


figure; hold on;
listchunk = floor(length(avocells) / 10);
lvl = 0;
for i = 1: listchunk:length(avocells)
  stop = i+(listchunk-1);
  if stop > length(avocells),
    stop = length(avocells);
  end
  plot(VPr_t, nanmean(VPr_Rpassive(avocells(i:stop),:))+lvl); hold on;
  lvl = lvl+1;
end
title('avocado')


figure; hold on;
listchunk = floor(length(lemcells) / 10);
lvl = 0;
for i = 1: listchunk:length(lemcells)
  stop = i+(listchunk-1);
  if stop > length(lemcells),
    stop = length(lemcells);
  end
  plot(VPr_t, nanmean(VPr_Rpassive(lemcells(i:stop),:))+lvl); hold on;
  lvl = lvl+1;
end
title('lemon')


figure; hold on;
lvl = 0;
for i = 1: length(wascells)
  plot(VPr_t, VPr_Rpassive(wascells(i),:)+lvl); hold on;
  lvl = lvl+2;
end
title('wasabi')

return;


