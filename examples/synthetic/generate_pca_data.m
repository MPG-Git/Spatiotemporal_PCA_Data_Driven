function generate_pca_data(destination)
if nargin<1, destination=fullfile(fileparts(mfilename('fullpath')),'data'); end
assert(~isfolder(destination),'Choose a new destination'); mkdir(destination);
oldRng=rng; cleanup=onCleanup(@()rng(oldRng)); rng(7);
one=fullfile(destination,'group_one'); two=fullfile(destination,'group_two'); mkdir(one); mkdir(two);
names={'F3','F4','C3','C4','P3','P4','Fz','Pz'};
locs=struct('labels',names,'theta',num2cell([-40 40 -80 80 -140 140 0 180]),'radius',num2cell([.35 .35 .4 .4 .35 .35 .2 .2]));
locs=convertlocs(locs,'topo2all');
for s=1:8
 E=eeg_emptyset; E.setname=sprintf('Synthetic subject %02d',s); E.subject=sprintf('synthetic_%02d',s); E.srate=100; E.nbchan=8; E.pnts=2000; E.trials=1;
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
 demo_portable_set(fullfile(folder,sprintf('s%02d.set',s)));
end

end
