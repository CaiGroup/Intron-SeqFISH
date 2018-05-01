

PathName = uigetdir;


folders = dir([PathName '/pos*']);

for i = 1:length(folders)
    num(i) = str2double(folders(i).name(4:end));
end

keep = num <30;

folders = folders(keep);

 %field = {'cpne5'; 'nes'; 'acta2';'gja1';'omg'; 'nov'; 'col5a1'; 'dcx';'itpr2';'rhob';'sox2';'cldn5'; 'mrc1'; 'tbr1';'pax6';...
    % 'calb1';'gda';'slc5a7'; 'sema3e'; 'mfge8';'lyve1';'loxl1';'slco1c1';'amigo2';'kcnip'};
    load([PathName '\' folders(1).name '\' folders(1).name 'Barcodesrad3ScaledOn.mat'])
    intron = copynumfinalrevised;
    positions = PosList;
    seedlist = seeds;
    Allint = int;
    Allints = ints;
    %cortotal = corr_offset;
     ba(1,1)=size(PosList,2)-1;
     ba(1,2)=size(copynumfinalrevised,2)-1;
     field = ones(1,size(copynumfinalrevised,2)-1)*str2double(folders(1).name(4:end));
      fullpath = [PathName '/' folders(1).name '/RoiSet'];
    vertex = selfseg(fullpath);
    for i = 1:length(vertex)
        BW = poly2mask(vertex(i).x, vertex(i).y, 2048, 2048);
        n = regionprops(BW,'Area','Centroid');
        area(i) = sum([n.Area]);
        centroid(i,:) = n.Centroid;
    end
    a = area;
    cent = centroid;
    clear area
    clear centroid
    
for i = 2:length(folders)
    load([PathName '\' folders(i).name '\' folders(i).name 'Barcodesrad3ScaledOn.mat'])
    %temp = copynumfinalrevised(:,2:end);
    intron = [intron, copynumfinalrevised(:,2:end)];
    positions = [positions, PosList(:,2:end)];
    seedlist = [seedlist seeds(:,2:end)];
    Allint = [Allint int(:,2:end)];
    Allints = [Allints ints(:,2:end)];
    field = [field,ones(1,size(copynumfinalrevised,2)-1)*str2double(folders(i).name(4:end))];
      fullpath = [PathName '/' folders(i).name '/RoiSet'];
    vertex = selfseg(fullpath);
    for j = 1:length(vertex)
        BW = poly2mask(vertex(j).x, vertex(j).y, 2048, 2048);
        n = regionprops(BW,'Area','Centroid');
        area(j) = sum([n.Area]);
        centroid(j,:) = n.Centroid;
    end
    cent = [cent;centroid];
    a = [a,area];
    clear area
    clear centroid
    %cortotal = [cortotal;corr_offset];
%     field.Positions = [field.Positions,temp.Positions];
%     field.Intensity = [field.Intensity,temp.Intensity];
%     temp = [];
     ba(i,1)=size(PosList,2)-1;
     ba(i,2)=size(copynumfinalrevised,2)-1;
end

