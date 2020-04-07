function beam=fcm_beam2Pimage(beam,pop,method,vFWHM,recalcs)
% FCM_BEAM2PIMAGE  creates P-images
%
% beam=fcm_beam2Pimage(beam_spatnorm,pop_controls,¦method¦,¦vFWHM¦,¦recalc_s¦)
%              (parameters in ¦¦ are optional)
%
% Input:
% beam_spatnorm     The spatially normalized (!) beam structure of a patient
% pop_controls      The values of all controls (pop structure)
% method            'normcdf' to calculate the position of the patient
%                     within the normal distribution of all controls at each 
%                     voxel (default).
%                   'ecdf' to calculate the position in the empirical
%                     distribution of all controls at each voxel.
% vFWHM             Gaussian kernel width in mm for variance smoothing.
%                   Default is [20 20 20]. Set to 0 or [] if you do not
%                   want to smooth the variance. Applies only to normcdf
%                   method.
% recalc_s          set this to true if you want to calculate a new ratio
%                   of the values in beam.s. By default, beam.s{1} is taken
%                   for calculating the P-image.
%
% Output:
% beam              beam structure with P-image data, which can be displayed 
%                   in nut_results_viewer.

if nargin<3 || isempty(method)
    method = 'normcdf';
end
if nargin<4 || isempty(vFWHM)
    vFWHM = [20 20 20];
elseif isscalar(vFWHM)
    vFWHM = vFWHM .* ones(1,3);
end
if ischar(beam)
    filename=beam;
    beam=load(filename);
    beam=nut_beam_legacy_compatibility(beam);
end
if ischar(pop)
    load(pop);
end
if length(beam.s)>1 && nargin>4 && recalcs
    [beam,selnum] = nut_calc_beamratio(beam);
else
    selnum=1;
end

numcontrolsubj = size(pop.s,1);
[beam.voxels,ia,ib] = intersect(beam.voxels,pop.voxels,'rows');  % find common voxels

switch lower(method)
    case 'normcdf'  

        % Get position in normal distribution of controls
        mu = shiftdim(mean(pop.s(:,ib,:,:),1),1);
        V = shiftdim(var(pop.s(:,ib,:,:),0,1),1);

        beam.s{1} = beam.s{1}(ia,:,:);   

        if ~isempty(vFWHM) && all(vFWHM)
            V = nut_smoothvol(V,beam,vFWHM);
        end 

        % p values from position in (normal) distribution of control subjects
        sigma = sqrt(V);
        
        if min(beam.s{1}(:))<0 && max(beam.s{1}(:))>0     
            beam.normcdf.p_uncorr = 2 * normcdf(-abs(beam.s{1}),mu,sigma);
            beam.normcdf.tail = 'both';
        elseif min(beam.s{1}(:))<0 && max(beam.s{1}(:))<1e-17
            beam.normcdf.p_uncorr = normcdf(beam.s{1},mu,sigma);
            beam.normcdf.tail = 'neg';
        elseif min(beam.s{1}(:))>-1e-17 && max(beam.s{1}(:))>0
            beam.normcdf.p_uncorr = normcdf(-beam.s{1},-mu,sigma);
            beam.normcdf.tail = 'pos';
        end
        beam.normcdf.p_uncorr(~isfinite(beam.normcdf.p_uncorr)) = 1;
        beam.normcdf.method='normcdf';
        
        % correct for multiple comparison (FDR)
%       now done in nut_FDR        
%         beam.normcdf.FDR=0.1;
%         [dum,beam.normcdf.cutoff]=nut_FDR(beam.normcdf.p_uncorr,[],beam.normcdf.FDR);
%         beam.normcdf.tail='both';    

        beam.s{1}(~isfinite(beam.s{1})) = 0;        % set eventual NaN's to 0   
        beam.s(2:end) = [];
        
        % Update labels
        if isfield(beam,'sinfo')
            beam.sinfo = beam.sinfo(selnum);
        end
        if isfield(beam,'labels') && isfield(beam.labels,'contrasts')
            beam.labels.contrasts = beam.labels.contrasts(selnum);
        end        
        
    case 'ecdf'
        beam.s{1} = beam.s{1}(ia,:,:); 
        I = shiftdim(beam.s{1},-1);
        P = pop.s(:,ib,:,:);
        
        if min(beam.s{1}(:))<0 && max(beam.s{1}(:))>0  
            N = sum(abs(P) > repmat(abs(I),[numcontrolsubj 1 1 1]),1);
            beam.normcdf.tail = 'both';
        elseif min(beam.s{1}(:))<0 && max(beam.s{1}(:))<1e-17
            N = sum(P < repmat(I,[numcontrolsubj 1 1 1]),1);
            beam.normcdf.tail = 'neg';
        elseif min(beam.s{1}(:))>-1e-17 && max(beam.s{1}(:))>0
            N = sum(P > repmat(I,[numcontrolsubj 1 1 1]),1);
            beam.normcdf.tail = 'pos';
        end
        N = shiftdim(N,1);
        
        beam.normcdf.p_uncorr = N./numcontrolsubj;
        beam.normcdf.p_uncorr(~isfinite(beam.normcdf.p_uncorr)) = 1;
        beam.normcdf.method='ecdf';
        
        beam.s{1}(~isfinite(beam.s{1})) = 0;        % set eventual NaN's to 0   
        beam.s(2:end) = [];
        
        if isfield(beam,'sinfo')
            beam.sinfo = beam.sinfo(selnum);
        end
        if isfield(beam,'labels') && isfield(beam.labels,'contrasts')
            beam.labels.contrasts = beam.labels.contrasts(selnum);
        end                
        
    case {'t' 'ttest' 't-test'}
        error('ttest method is strongly discouraged because of false positives.')
%         % T-tests, 0 hypothesis: difference between patient and controls is 0.
%         temp(1,:,:,:) = beam.s{1}(ia,:,:);
%         test = repmat(temp,[numcontrolsubj 1 1 1]) - pop.s(:,ib,:,:);
%         V = var(test,0,1);
%         if ~isempty(vFWHM) && all(vFWHM)
%             V(1,:,:,:) = nut_smoothvol(shiftdim(V,1),beam,vFWHM);
%         end 
%         
%         T = shiftdim( mean(test,1) ./ sqrt(  V ./ numcontrolsubj ) ,1 );
%         p = 2.*tcdf(-abs(T), numcontrolsubj-1);
% 
%         T(~isfinite(T)) = 0;        % set eventual NaN's to 0
%         p(~isfinite(p)) = 1;
% 
%         % The "power change" values displayed in nut_timef_viewer are T values
%         beam.s{1}=T;
%         beam.sinfo={'T'};
%         beam.s(2:end)=[];
% 
%         % p values from ttest
%         beam.ttest.T = T;
%         beam.ttest.tail='both';  
%         beam.ttest.p_uncorr = p;
%         beam.ttest.FDR=0.01;
%         
%         % Update labels
%         beam.sinfo = {'T'};
%         if isfield(beam,'labels')
%             if isfield(beam.labels,'contrasts')
%                 beam.labels = rmfield(beam.labels,'contrasts');
%             end
%             if size(beam.s{1},3)>1
%                 beam.labels.colorbar = 'T';
%             else
%                 beam.labels.yaxis = 'T';
%             end       
%         end
%         
    case 'snpm'
        % T-tests, 0 hypothesis: difference between patient and controls is 0.
        temp(1,:,:,:) = beam.s{1}(ia,:,:);
        pop.s = repmat(temp,[numcontrolsubj 1 1 1]) - pop.s(:,ib,:,:);
        pop.voxels=pop.voxels(ib,:);
        save popTemp pop
        
        stats = struct('type','ttest1_snpm','docomp',false,'fsel','one','frqstats',false, ...
            'corr4multfreq',false,'tsel','one','tsstats',false,'corr4multtime',false, ...
            'covar',[],'whichtail','both','markmode',2,'doavg',true,'dostat',true, ...
            'usespatnormed',true,'subjsel',unique(pop.subjnr),'numsubj',length(unique(pop.subjnr)), ...
            'condsel',unique(pop.condnr),'numtests',1,'comps',1,'subjnr',pop.subjnr, ...
            'condnr',pop.condnr,'files',{{'popTemp.mat'}}, ...
            'frqband',[],'frqtxt','','timepts',[],'timetxt','');
        stats.numperm =  2^stats.numsubj;
        stats.isaproxperm = (stats.numperm>1e4);
        if stats.isaproxperm, stats.numperm=1e4; end
        
        nut_beamstats(stats);
        delete popTemp.mat 
        beam=load('s_beamtf1_ttest1');
%         delete s_beamtf1_avg.mat s_beamtf1_ttest1.mat
        
    otherwise
        error('Unknown method for creating P-images.')
end

 % update MRI coreg info
if strncmp( spm('ver'),'SPM8',4 ) 
    spmmri = [fileparts(which('spm')) filesep 'canonical' filesep 'avg152T1.nii'];
elseif strcmp( spm('ver'),'SPM2' )
    spmmri = [fileparts(which('nutmeg')) filesep 'templates' filesep 'wNormT1neuroax.img'];
else
    error('This version of SPM is currently not supported.')
end
beam.coreg = struct('mripath', spmmri, ...
     'meg2mri_tfm',eye(4),'orientation',1, ...
     'norm_mripath',spmmri, ...
     'brainrender_path',[fileparts(which('spm.m')) filesep 'rend' filesep 'render_single_subj.mat']);

 if nargout<1 && exist('filename','var')
     [pa,fi,ext]=fileparts(filename);
     fi=strrep(fi,'_spatnorm','P');
     save(fullfile(pa,fi),'beam');
 end



