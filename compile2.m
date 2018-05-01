
PathName = uigetdir;


folders = dir([PathName '/pos*']);

for i = 1:length(folders)
    num(i) = str2double(folders(i).name(4:end));
end

keep = num <30;

folders = folders(keep);

load([PathName '\' folders(1).name '\Nucleolus.mat'])
nuc = nucleolus;
    
for i = 2:length(folders)
    load([PathName '\' folders(i).name '\Nucleolus.mat'])
    nuc = [nuc, nucleolus];
end

