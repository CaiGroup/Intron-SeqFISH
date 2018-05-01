load('G:\Yodai\left_confocal\20180117_Yodai_7_short_E14_serum_10kintron_rep3_trial2\decodeVars.mat');
for i = 0:14 
    load(['G:\Yodai\right_confocal\20180108_NIH3T3_rep1\Organized\Pos' num2str(i) '\pos' num2str(i) 'PointsBarcodesrad6ScaledOn.mat']);
    load(['G:\Yodai\right_confocal\20180108_NIH3T3_rep1\Organized\Pos' num2str(i) '\pos' num2str(i) 'Barcodesrad6ScaledOn.mat']);
    [dotlocations,copynumfinalrevised,PosList] = PointLocations('G:\Yodai\right_confocal\20180108_NIH3T3_rep1\Organized',i,5, 12, points, foundbarcodes,barcodekey,copynumfinal,3); 
    [seeds, int, ints] = numseeds(PosList,dotlocations,5);
    save(['G:\Yodai\right_confocal\20180108_NIH3T3_rep1\Organized\Pos' num2str(i) '\pos' num2str(i) 'SeedsInt.mat'],'dotlocations','seeds','ints','int'); 
end