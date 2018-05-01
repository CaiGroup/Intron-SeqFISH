function save_for_parfor(PathName,field,foundbarcodes,points,copynumfinal,copynumfinalsum,totdropped,rawfound,PosList,dotlocations,copynumfinalrevised, seeds, int, ints)

save([PathName '\' 'pos' num2str(field) '\' 'pos' num2str(field) 'Barcodesrad3ScaledOn.mat'],'foundbarcodes', 'copynumfinal', 'copynumfinalsum', 'totdropped','rawfound','PosList','dotlocations','copynumfinalrevised','seeds','int','ints','-v7.3');

save([PathName '\' 'pos' num2str(field) '\' 'pos' num2str(field) 'PointsBarcodesrad3ScaledOn.mat'],'points','-v7.3')
    clear foundbarcodes
    clear points
    clear copynumfinal
    clear rawfound
    clear hybnum
    clear copynumfinalsum
    clear totdropped
    clear dotlocations
    clear copynumfinalrevised
    clear PosList
    clear int
    clear seeds
    clear ints