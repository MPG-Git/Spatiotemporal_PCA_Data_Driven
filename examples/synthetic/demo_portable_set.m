function demo_portable_set(file)
% Keep generated .set files independent of the author's directory layout.
stored=load(file,'-mat'); stored.filepath='';
save(file,'-struct','stored','-mat');
end
