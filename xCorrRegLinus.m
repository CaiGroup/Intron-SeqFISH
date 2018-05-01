function [hybnum, corr_offset, images] = xCorrRegLinus( PathName,hybnum,hyb,channels,posnum )
fld = pwd;
Miji;
cd(fld);
%find offsets
for i = 1:hyb
    path_to_fish = ['path=[' PathName '\' 'pos' num2str(posnum) '\hyb' num2str(i) 'RegistrationCheck.tif' ']'];

    MIJ.run('Open...', path_to_fish);
    MIJ.run('Split Channels')
    if i == 1
        images{1} = MIJ.getImage('C1-hyb1RegistrationCheck.tif');
        MIJ.run('Close All')
        if size(images{1},3)<16
            add = repmat(images{1}(:,:,end),1,1,16-size(images{1},3));
            images{1} = cat(3,images{1},add);
        end
    else
        imagestemp = MIJ.getImage(['C1-hyb' num2str(i) 'RegistrationCheck.tif']);
        MIJ.run('Close All')
        if size(imagestemp,3)<16
            add = repmat(imagestemp(:,:,end),1,1,16-size(imagestemp,3));
            imagestemp = cat(3,imagestemp,add);
        end
        [optimizer, metric] = imregconfig('monomodal');
        tform = imregtform(imagestemp,images{1},'translation',optimizer,metric);
        images{i} = imwarp(imagestemp,tform,'OutputView',imref3d(size(imagestemp)));
        corr_offset{i} = tform;
        for j = 1:length(channels)
            hybnum(i).color{j} = imwarp(hybnum(i).color{j},tform,'OutputView',imref3d(size(hybnum(i).color{j})));
        end
    end
end

for i = 1:hyb
    zend(i) = size(images{i},3);
end

stop = min(zend);

for i = 1:hyb
    MIJ.createImage(num2str(i),images{i}(:,:,1:stop),true)
end
str = [];
for i = 1:hyb
    str = [str ' image' num2str(i) '=' num2str(i)];
end
MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);      
MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(hyb) ' slices=' num2str(stop) ' frames=1 display=Grayscale']);
MIJ.run('Save', ['save=[' PathName '\Pos' num2str(posnum) '\AllHybRegistrationCheck.tif' ']']);
MIJ.run('Close All')
MIJ.exit

