function smoke_spatiotemporal_pca
addpath(fileparts(fileparts(mfilename('fullpath'))));
assert(~isempty(which('eeg_emptyset')),'Initialize EEGLAB first');
set(groot,'defaultFigureVisible','off'); rng(7);
root=tempname; mkdir(root);
one=fullfile(root,'one'); two=fullfile(root,'two');
if ~isfolder(one), mkdir(one); mkdir(two); end
names={'F3','F4','C3','C4','P3','P4','Fz','Pz'};
locs=struct('labels',names,'theta',num2cell([-40 40 -80 80 -140 140 0 180]),'radius',num2cell([.35 .35 .4 .4 .35 .35 .2 .2]));
locs=convertlocs(locs,'topo2all');
for s=1:8
 E=eeg_emptyset; E.srate=100; E.nbchan=8; E.pnts=2000; E.trials=1;
 E.xmin=0; E.xmax=19.99; E.times=(0:1999)*10; E.chanlocs=locs;
 E.data=randn(8,2000)*.2; amps=5+randn(8,1)*2; amps2=randn(8,1)*5;
 lat=150:250:1800; E.event=struct('type',{},'latency',{},'duration',{});
 for k=1:numel(lat)
  E.data=E.data+amps*exp(-(((1:2000)-lat(k)-30)/6).^2);
  E.data=E.data+amps2*exp(-(((1:2000)-lat(k)-50)/4).^2);
  E.event(k)=struct('type','stim','latency',lat(k),'duration',0);
 end
 E=eeg_checkset(E,'eventconsistency');
 if s<=4, folder=one; else, folder=two; end
 pop_saveset(E,'filename',sprintf('s%02d.set',s),'filepath',folder,'savemode','onefile');
end
E.chanlocs=E.chanlocs([2 1 3:8]);
pop_saveset(E,'filename','z_mismatched.set','filepath',two,'savemode','onefile');
before=path;
run_spatiotemporal_pca_wrapper(one,two,fullfile(root,'results'),{'stim'},[-.2 .7],[-200 0],[200 400],10,{{'varimax'},{'promax','Power',4}});
assert(strcmp(before,path),'The MATLAB path changed');
files=dir(fullfile(root,'results','pca_*.mat')); assert(~isempty(files));
assert(contains(fileread(fullfile(root,'results','mismatch_length_all.txt')),'z_mismatched.set'));
for k=1:numel(files)
 result=load(fullfile(files(k).folder,files(k).name));
 original=result.score_time_unrot*result.coeff_time_unrot';
 rotated=result.rotated_scores*result.coeff_time_rot';
 assert(norm(original-rotated,'fro')/norm(original,'fro')<1e-10);
end
assert(~isempty(dir(fullfile(root,'results','*.png'))));
fprintf('PCA SYNTHETIC DATA, ROTATIONS, RECONSTRUCTION, OUTPUTS, PATH: PASS\n');
close all;
end
