function replaceChannel(PathName1, PathName2, channel1, channel2, replace, imchan, z)
listing1 = dir([PathName1 '\*.tif']);
listing2 = dir([PathName2 '\*.tif']);

fld = pwd;
Miji;
cd(fld);

for i = 1:length(listing1)
    
    path_to_fish = ['path=[' PathName1 '\' listing1(i).name ']'];
    MIJ.run('Open...', path_to_fish);
    MIJ.run('Split Channels')
    for channum = 1:length(channel1)
        im{channum} = uint16(MIJ.getImage(['C' num2str(channel1(channum)) '-' listing1(i).name]));
    end
    MIJ.run('Close All')
    
    path_to_fish = ['path=[' PathName2 '\' listing2(i).name ']'];
    MIJ.run('Open...', path_to_fish);
    MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(imchan) ' slices=' num2str(z) ' frames=1 display=Grayscale']);
    MIJ.run('Split Channels')
    for channum = 1:length(channel2)
        im2{channum} = uint16(MIJ.getImage(['C' num2str(channel2(channum)) '-' listing2(i).name]));
    end
    MIJ.run('Close All')
    
%     tform = imregcorr(max(im2{end},[],3),max(im{end},[],3));
%     for channum = 1:length(im2)
%         im2{channum} = imwarp(im2{channum},tform,'OutputView',imref2d(size(im{end}(:,:,1))));
%     end
    
    [optimizer, metric] = imregconfig('monomodal');
    tform = imregtform(im2{end},im{end},'translation',optimizer,metric);
%     tform.T(4,1) = 0;
%     tform.T(4,2) = 0;
    for channum = 1:length(im2)
        im2{channum} = imwarp(im2{channum},tform,'OutputView',imref3d(size(im{end})));
    end
    
    if size(im2{1},3)>size(im{1},3)
        for channum = 1:length(im2)
            im2{channum} = im2{channum}(1:size(im{1},3));
        end
    elseif size(im2{1},3)>size(im{1},3)
        for channum = 1:length(im1)
            im1{channum} = im1{channum}(1:size(im2{1},3));
        end
    end
    
    MIJ.createImage('1', im{end}, true);
    MIJ.createImage('2', im2{end}, true);
    
    MIJ.run('Merge Channels...', 'c1=1 c2=2 create');
    
    MIJ.run('Save', ['save=[' PathName1 '\Reg' listing1(i).name ']']);
    
    MIJ.run('Close All')
    
    for channum = 1:length(channel1)
        MIJ.createImage(['Orig' num2str(channel1(channum))], im{channum},true)
    end
    
    for channum = 1:length(channel2)-1
        MIJ.createImage(['replace' num2str(replace(channum))], im2{channum},true)
    end
    
    str = [];
    for channum = 1:4
        if ismember(channum,replace)
            str = [str ' c' num2str(channum) '=replace' num2str(replace)];
        else
            str = [str ' c' num2str(channum) '=Orig' num2str(channum)];
        end
    end
    
    MIJ.run('Merge Channels...', [str ' create']);
    
    MIJ.run('Save', ['save=[' PathName1 '\' listing1(i).name ']']);
    
    MIJ.run('Close All')
    
end
MIJ.exit
