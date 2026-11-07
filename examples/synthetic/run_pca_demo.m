function run_pca_demo(out)
if nargin<1, out=''; end
[here,out,cleanup]=demo_prepare(mfilename('fullpath'),out); %#ok<ASGLU>
addpath(fileparts(fileparts(here))); rng(71);
before=findall(groot,'Type','figure');
run_spatiotemporal_pca_wrapper(fullfile(here,'data','group_one'), ...
 fullfile(here,'data','group_two'),out,{'stim'},[-.2 .7],[-200 0],[200 400], ...
 10,{{'varimax'},{'promax','Power',4}});
created=setdiff(findall(groot,'Type','figure'),before);
fc=onCleanup(@()close(created(isgraphics(created)))); %#ok<NASGU>
files=dir(fullfile(out,'pca_*.mat')); assert(~isempty(files));
names={}; values=[];
for k=1:numel(files)
 r=load(fullfile(files(k).folder,files(k).name));
 original=r.score_time_unrot*r.coeff_time_unrot'; rotated=r.rotated_scores*r.coeff_time_rot';
 err=norm(original-rotated,'fro')/norm(original,'fro'); assert(err<1e-10);
 names=[names,{[files(k).name ':TemporalFactors'],[files(k).name ':SpatialFactors'],[files(k).name ':ReconstructionError']}]; %#ok<AGROW>
 values=[values,size(r.coeff_time_rot,2),size(r.coeff_space_rot,2),err]; %#ok<AGROW>
 % Sign-invariant, numeric output for independent inspection.
 gram=r.coeff_time_rot*r.coeff_time_rot';
 gramName=[files(k).name '_temporal_gram.csv'];
 writematrix(gram,fullfile(out,gramName));
 reference=readmatrix(fullfile(here,'expected',gramName));
 assert(isequal(size(gram),size(reference)) && ...
  norm(gram-reference,'fro')/max(norm(reference,'fro'),eps)<1e-6, ...
  'Temporal loading Gram matrix differs from the bundled reference');
end
assert(~isempty(dir(fullfile(out,'*.png'))));
demo_verify(here,out,names,values);
end
