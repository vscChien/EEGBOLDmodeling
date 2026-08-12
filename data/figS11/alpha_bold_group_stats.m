function alpha_bold_group_stats()
% ALPHA_BOLD_GROUP_STATS  Group statistics and FDR-thresholded map of the
% empirical alpha-BOLD correlation (Chien et al., Fig S11).
%
% Two-tailed one-sample sign-flip permutation test of the mean correlation
% against zero, Benjamini-Hochberg FDR correction across voxels, and a
% BrainNet Viewer rendering of the surviving voxels as PNG.
%
% USAGE
%   Add FieldTrip and BrainNet Viewer to the MATLAB path, then run:
%       alpha_bold_group_stats
%
% DEPENDENCIES (not bundled - add to path yourself)
%   FieldTrip           https://www.fieldtriptoolbox.org
%                       used: ft_defaults, ft_sourceinterpolate, ft_volumewrite
%   BrainNet Viewer     https://www.nitrc.org/projects/bnv
%                       used: BrainNet_MapCfg
%
% BUNDLED DATA (in ./data, loaded by relative path)
%   AlphaBOLDCorr.mat   AlphaBOLDCorr [4417 x 72] Pearson r per voxel/subject
%                       inside, pos, unit - source grid definition
%   mni_sourcemodel.mat 6 mm grey matter source grid (cm)
%   mni_mri.mat         MNI template used as interpolation target
%   mycmap.mat          colour map for BrainNet
%   BNopt1.mat          BrainNet option struct
%
% OUTPUT (written next to this script)
%   alpha_bold_fdr.nii            thresholded map as NIfTI
%   alpha_bold_fdr_q0p05.png      rendered figure
%   alpha_bold_group_stats.mat    meanR, p, q, mask
%
% Runtime ~2 min for 1e5 permutations on 4417 voxels.

%% Settings
nPerm = 100000;   % sign-flip permutations
qFDR  = 0.05;     % false discovery rate level
seed  = 42;       % for reproducibility

%% Paths - everything relative to this file
here = fileparts(mfilename('fullpath'));
dataDir = fullfile(here,'data');

assert(exist('ft_defaults','file')==2, ...
    'FieldTrip not on the MATLAB path.')
assert(exist('BrainNet_MapCfg','file')==2, ...
    'BrainNet Viewer not on the MATLAB path.')
ft_defaults

load(fullfile(dataDir,'AlphaBOLDCorr.mat'),'AlphaBOLDCorr');
load(fullfile(dataDir,'mni_sourcemodel.mat'),'mni_sourcemodel');
load(fullfile(dataDir,'mni_mri.mat'),'mni_mri');
load(fullfile(dataDir,'mycmap.mat'),'mycmap');

X = AlphaBOLDCorr;                 % [nVox x nSub] raw Pearson r
X(isnan(X)) = 0;
[nVox,nSub] = size(X);
meanR = mean(X,2);

fprintf('%d voxels, %d subjects, %d permutations\n',nVox,nSub,nPerm);

%% Sign-flip permutation null
% One matrix product gives a block of permutations for all voxels at once.
rng(seed)
S = sign(randn(nSub,nPerm));       % observed labelling deliberately not in S
cnt = zeros(nVox,1);
absObs = abs(meanR);

blk = 500;
for a = 1:blk:nPerm
    b   = min(a+blk-1,nPerm);
    muP = (X * S(:,a:b))/nSub;
    cnt = cnt + sum(abs(muP) >= absObs, 2);
end

% Two-tailed p, +1 correction so that p is never zero (Phipson & Smyth 2010)
p = (cnt + 1)/(nPerm + 1);

%% Benjamini-Hochberg FDR
[mask,qval] = bh_fdr(p,qFDR);

fprintf('FDR q<%.2f: %d/%d voxels significant (%d positive, %d negative)\n', ...
    qFDR,sum(mask),nVox,sum(mask & meanR>0),sum(mask & meanR<0));

save(fullfile(here,'alpha_bold_group_stats.mat'), ...
    'meanR','p','qval','mask','nPerm','qFDR','seed');

%% Map the thresholded values onto the full grid and write NIfTI
data = meanR .* double(mask);      % non-significant voxels -> 0

grid_all = mni_sourcemodel;
grid_all.inside(:) = true;
d = pdist2(mni_sourcemodel.pos(mni_sourcemodel.inside,:),grid_all.pos,'euclidean');
[~,nearest] = min(d,[],1);
grid_all.map = data(nearest)';

cfg = [];
cfg.parameter    = 'map';
cfg.interpmethod = 'nearest';
cfg.keepinside   = 'no';
vol = ft_sourceinterpolate(cfg,grid_all,mni_mri);

niiFile = fullfile(here,'alpha_bold_fdr');
cfg = [];
cfg.filename  = niiFile;
cfg.parameter = 'map';
cfg.filetype  = 'nifti';
cfg.datatype  = 'float';
ft_volumewrite(cfg,vol);

%% Render with BrainNet Viewer
load(fullfile(dataDir,'BNopt1.mat'),'EC');
EC.vol.color_map = 25;
EC.vol.CM_annot  = mycmap;
EC.vol.cmstring  = 'Custom';
EC.vol.px = max(data);             % positive extreme
EC.vol.nx = min(data);             % negative extreme
EC.vol.pn =  0;                    % no dead zone, the FDR mask defines
EC.vol.nn = -0;                    % what is shown
optFile = fullfile(here,'BNopt_fdr.mat');
save(optFile,'EC');

pngFile = fullfile(here,sprintf('alpha_bold_fdr_q%s.png', ...
    strrep(sprintf('%g',qFDR),'.','p')));
BrainNet_MapCfg('BrainMesh_ICBM152_tal.nv',[niiFile '.nii'],optFile);
print(pngFile,'-dpng','-r300')
close all

fprintf('wrote %s\n',pngFile);
end

%% ------------------------------------------------------------------------
function [h,qval] = bh_fdr(p,q)
% Benjamini-Hochberg step-up FDR. Returns the rejection mask and the
% adjusted p-values (q-values).
p = p(:);
m = numel(p);
[ps,idx] = sort(p);
thresh = (1:m)'/m * q;
last = find(ps <= thresh, 1, 'last');

h = false(m,1);
if ~isempty(last)
    h(idx(1:last)) = true;
end

% adjusted p-values: running minimum from the largest p downwards
adj = min(1, cummin((m./(m:-1:1)') .* flipud(ps)));
qval = zeros(m,1);
qval(idx) = flipud(adj);
end
