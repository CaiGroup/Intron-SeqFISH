function offsets = make12channel3D(PathName, FolderName, hybs, channelstotal, channelsperhyb)

for i = 1:hybs
    listing{i} = dir([PathName '\' num2str(i) '\*.tif']);
end
mkdir(PathName,FolderName);
fld = pwd;
Miji;
cd(fld);
k2 = strfind(listing{1}(1).name, 'Pos');
        
for k = 1:length(listing{1})
    roi = listing{1}(k).name(k2(1):k2(1)+4);
    mkdir([PathName '\' FolderName],roi);
    for i = 0:(hybs/(channelstotal/channelsperhyb))-1
        counter = 1;
        for j = (i*channelstotal/channelsperhyb)+1:channelstotal/channelsperhyb*(i+1)
            Path = ['path=[' PathName '\' num2str(j) '\' listing{j}(k).name ']'];
            MIJ.run('Open...', Path);
            %hybtemp = uint16(MIJ.getCurrentImage);
            MIJ.run('Split Channels')
            for channum = 1:channelsperhyb+1
                im{channum} = uint16(MIJ.getImage(['C' num2str(channum) '-' listing{j}(k).name]));
            end
            MIJ.run('Close All')
            if mod(j,channelstotal/channelsperhyb) ~=1
                regi = cat(2,regi,im{4});
%                 c = normxcorr2(regi(:,:,counter),regi(:,:,1));
%                 [~, imax] = max(abs(c(:)));
%                 [ypeak, xpeak] = ind2sub(size(c),imax(1));
%                 offset{i+1,counter} = [(xpeak-size(regi(:,:,1),2)) (ypeak-size(regi(:,:,1),1))];
%                 hybtempreg = imtranslate(hybtemp,offset{i+1,counter});
                [optimizer, metric] = imregconfig('monomodal');
                tform = imregtform(regi{counter},regi{1},'translation',optimizer,metric);
                for channum = 1:channelsperhyb+1
                    imreg{channum} = imwarp(im{channum},tform,'OutputView',imref3d(size(regi{1})));
                end
                hyb = cat(2,hyb, imreg(1:3));
                regireg = cat(2,regireg, imreg{4});
                offsets{i+1,counter} = tform;
            else
                hyb = im(1:3);
                regireg = im(4);
                regi = im(4);
                hybnum(1).color{1} = im{4};
                %zstart = AutoFocus(hybnum, 1, 1,0)-3;
                zstart = 1;
            end
            counter = counter + 1;
            MIJ.run('Close All')
        end
        for channum = 1:channelstotal 
            MIJ.createImage(num2str(channum),hyb{channum}(:,:,zstart:end),true); 
        end
        str = [];
        for channum = 1:channelstotal
            str = [str ' image' num2str(channum) '=' num2str(channum)];
        end
        MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);      
        MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(channelstotal) ' slices=' num2str(size(hyb{1}(:,:,zstart:end),3)) ' frames=1 display=Grayscale']);
        MIJ.run('Save', ['save=[' PathName '\' FolderName '\' roi '\' num2str(i+1) '.tif' ']']);
        MIJ.run('Close All')
        for channum = 1:channelstotal/channelsperhyb 
            MIJ.createImage(num2str(channum),regireg{channum}(:,:,zstart:end),true); 
        end
        str = [];
        for channum = 1:channelstotal/channelsperhyb
            str = [str ' image' num2str(channum) '=' num2str(channum)];
        end
        MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);      
        MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(channelstotal/channelsperhyb) ' slices=' num2str(size(hyb{1}(:,:,zstart:end),3)) ' frames=1 display=Grayscale']);
        MIJ.run('Save', ['save=[' PathName '\' FolderName '\' roi '\hyb' num2str(i+1) 'RegistrationCheck.tif' ']']);
        MIJ.run('Close All')
    end
end
        
MIJ.exit