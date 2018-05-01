function HyperSwap(PathName,zslice)
%PathName = 'G:\New Intron\20171216_Yodai_6_short_E14_serum_10kintron_2';
listing = dir([PathName '\*']);
listing(1:2) = [];
fld = pwd;
Miji;
cd(fld);
for i = 1:length(listing)
    listing2 = dir([PathName '\' listing(i).name '\*.tif']);
    for j = 1:length(listing2)
        path_to_fish = ['path=[' PathName '\' listing(i).name '\' listing2(j).name ']'];
        MIJ.run('Open...', path_to_fish);
        MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=4 slices=' zslice ' frames=1 display=Grayscale']);
        MIJ.run('Arrange Channels...', 'new=1324');
        MIJ.run('Save', ['save=[' PathName '\' listing(i).name '\' listing2(j).name ']']);
        MIJ.run('Close All')
    end
end
MIJ.exit