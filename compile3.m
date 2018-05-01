PathName = uigetdir;


folders = dir([PathName '/pos*']);

for i = 1:length(folders)
    num(i) = str2double(folders(i).name(4:end));
end

keep = num <30;

folders = folders(keep);

 %field = {'cpne5'; 'nes'; 'acta2';'gja1';'omg'; 'nov'; 'col5a1'; 'dcx';'itpr2';'rhob';'sox2';'cldn5'; 'mrc1'; 'tbr1';'pax6';...
    % 'calb1';'gda';'slc5a7'; 'sema3e'; 'mfge8';'lyve1';'loxl1';'slco1c1';'amigo2';'kcnip'};
    load([PathName '\' folders(1).name '\' folders(1).name 'SeedsInt6.mat'])
    seedlist = seeds;
    Allint = int;
    Allints = ints;
    PL = PosList;
    cpr = copynumfinalrevised;
    
for i = 2:length(folders)
    load([PathName '\' folders(i).name '\' folders(i).name 'SeedsInt6.mat'])
    seedlist = [seedlist seeds(:,2:end)];
    Allint = [Allint int(:,2:end)];
    Allints = [Allints ints(:,2:end)];
    PL = [PL PosList(:,2:end)];
    cpr = [cpr copynumfinalrevised(:,2:end)];
end

