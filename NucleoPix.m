function [nucleolus]= NucleoPix(PathName, imnum, channel, regchan, regist)

fld = pwd;
Miji;
cd(fld);

path_to_fish = ['path=[' PathName '\' num2str(imnum) '.tif]'];
MIJ.run('Open...', path_to_fish);
MIJ.run('Split Channels');

its1 = uint16(MIJ.getImage(['C' num2str(channel) '-' num2str(imnum) '.tif']));
reg = uint16(MIJ.getImage(['C' num2str(regchan) '-' num2str(imnum) '.tif']));

MIJ.run('Close All')

path_to_fish = ['path=[' PathName '\' regist '.tif]'];
MIJ.run('Open...', path_to_fish);

fixed = uint16(MIJ.getCurrentImage);

MIJ.run('Close All')

MIJ.createImage('its1',its1,true)

MIJ.run('Auto Threshold', 'method=Default white stack use_stack_histogram');

its1 = uint16(MIJ.getCurrentImage) >0;

MIJ.run('Close All')
MIJ.exit

[optimizer, metric] = imregconfig('monomodal');
tform = imregtform(reg,fixed,'translation',optimizer,metric);

its1 = imwarp(its1,tform,'OutputView',imref3d(size(fixed)));

nucleomask = bwperim(its1);
fullpath = [PathName '\RoiSet'];
vertex = selfseg(fullpath);

[y,x,z] = ind2sub(size(nucleomask),find(nucleomask == 1));

for i = 1:length(vertex)
    include = inpolygon(x,y,vertex(i).x,vertex(i).y);
    nucleolus{i} = [x(include) y(include) z(include)];
end



