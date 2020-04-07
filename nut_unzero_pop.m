function pop=nut_unzero_pop(pop,cut,idx3dim)
% pop=nut_unzero_pop(pop,cut,idx3dim)

if nargin<2 || isempty(cut), cut=0; end
if nargin<3, idx3dim=1; end

bad=(sum(pop.s(:,:,idx3dim)==0,1)>cut)
pop.s(:,bad,:)=[];
pop.voxels(bad,:)=[];