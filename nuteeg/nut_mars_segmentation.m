function nut_mars_segmentation
% NUT_MARS_SEGMENTATION calls external toolboxes to segment scalp, skull and brain.


global nuts

if exist('start_mars.m','file')<2 || exist('mysegment.m','file')<2 ||  exist('load_untouch_nii.m','file')<2
    errordlg(sprintf(['This needs the following toolboxes in the Matlab path: \n' ... 
         '1. MARS: https://www.parralab.org/mars/ \n' ...
         '2. Segment: http://neuralengr.com/segment/ \n' ...
         '3. ANALYZE: http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image. \n' ...
         'Also make sure you have the lastest update of SPM8 (use spm_update).']))
    return
end
    
if isempty(nuts)
    [file,path] = uigetfile('*.nii;*.img','Open MRI for segmentation');
    if isequal(file,0); return; end;
    mripath = fullfile(path,file);
else
    mripath = nuts.coreg.mripath;
end
[path,file,ext]=fileparts(mripath);


% start_mars(T1,T2,TPM,TCM,numGaussian,runningMode,beta,convergence,saveBiasField,saveNative,saveWarp,saveDeformField,biasRegul,biasFWHM)
%           The two last input variables are only on a modified start_mars.m file (ask Lea)
% Default: start_mars([],[],'eTPM.nii','rTCM_BW20_S1.mat',[2 2 2 3 4 2],'MARS',0.1,0.1,'Save Nothing','Native Space','None','None', 0.0001,60)
% Default: start_mars(mripath,[],[],[], [],[],0.1, 0.1,[],[], [],[], 0.0001,60);

start_mars(mripath);
mysegment(path,file);

L={'mask_gray' 'mask_white' 'mask_csf' 'mask_bone' 'mask_skin'};
l   = {'gray', 'white', 'csf','skull','scalp'};
mri = ft_read_mri(mripath);
seg = mri;
for k=1:5
    mri = ft_read_mri([L{k} '.nii']);
    seg.(l{k})=mri.anatomy;
end
clear mri

scalp = seg.scalp;
skull = seg.skull;
brain = seg.gray + seg.white + seg.csf;
spm_smooth(brain,brain,5);
brain = volumethreshold(brain, 0.5);

% To avoid tiny brain stem on caudal end, which causes errors in surface
numvoxz=shiftdim(sum(sum(brain,1),2));
beginhere=find(numvoxz>100,1);
brain(:,:,1:beginhere-1)=false;

scalp=logical(scalp); skull=logical(skull); 
for k=1:size(scalp,3)
    skull(:,:,k)=bwconvhull(skull(:,:,k));
    scalp(:,:,k)=bwconvhull(scalp(:,:,k));
end
scalp=double(scalp); skull=double(skull);

% the mex files underneath the spm function will change the input variable
spm_smooth(scalp, scalp, 20);
scalp = volumethreshold(scalp, 0.5);

spm_smooth(skull,skull,20);
skull = volumethreshold(skull, 0.5);

scalp(brain>0)=0;
scalp(skull>0)=0;
skull(brain>0)=0;

seg.scalp = scalp;
seg.skull = skull;   
seg.brain = brain;
seg = removefields(seg, {'anatomy' 'gray' 'white' 'csf'});

segi = ft_datatype_segmentation(seg,'segmentationstyle','indexed');
% for some reason, the order of the segments is sometime reversed.
% It seems MATLAB changed the behavior of intersect used in 
% ft_datatype_segmentation, which inverses the order of the segments!
% WTF MATLAB!
% We have to correct for this possibility here:
order={'scalp', 'skull', 'brain'};
if ~isequal(segi.seglabel,order)
    stmp=segi.seg;
    for n=1:3
        segi.seg(stmp==find(strcmp(segi.seglabel,order{n})))=n;
    end
    segi.seglabel=order;
end

ft_write_mri(fullfile(path,[file '.spmlabel.nii']),segi.seg,'dataformat','nifti','transform',segi.transform);

nut_iso2surf(fullfile(path,[file '.spmlabel.nii']));
% cfg = [];
% cfg.method='projectmesh';
% cfg.numvertices = [1000 2000 3000];
% cfg.tissue      = {'scalp', 'skull', 'brain'};
% bnd = ft_prepare_mesh(cfg, seg);



%----------------------------------------------
function [output] = volumethreshold(input, thresh)

% private function taken from FieldTrip.
%
% VOLUMETHRESHOLD is a helper function for segmentations. It applies a
% relative threshold and subsequently looks for the largest connected part,
% thereby removing small blobs such as vitamine E capsules.

% mask by taking the negative of the segmentation, thus ensuring
% that no holes are within the compartment and do a two-pass
% approach to eliminate potential vitamin E capsules etc.

if ~islogical(input)
  output = double(input>(thresh*max(input(:))));
else
  % there is no reason to apply a threshold, but spm_bwlabel still needs a
  % double input for clustering
  output = double(input);
end

% cluster the connected tissue
[cluster, n] = spm_bwlabel(output, 6);

if n>1
  % it pays off to sort the cluster assignment if there are many clusters
  tmp = cluster(:);                       % convert to a vector
  tmp = tmp(tmp>0);                       % remove the zeros
  tmp = sort(tmp, 'ascend');              % sort according to cluster number
  m   = zeros(1,n);
  for k=1:n
    m(k) = sum(tmp==k);       % determine the voxel count for each cluster
    tmp  = tmp(m(k)+1:end);   % remove the last cluster that was counted
  end
  % select the tissue that has the most voxels belonging to it
  [m, i] = max(m);
  output = (cluster==i);
else
  % the output only contains a single cluster
  output = (cluster==1);
end



