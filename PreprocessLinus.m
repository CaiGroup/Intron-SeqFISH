function [ hybnum] = PreprocessLinus( channels, hyb, PathName, posnum,tforms,corrections)
%Preprocess Summary of this function goes here
%   Detailed explanation goes here

fld = pwd;
Miji;
cd(fld);


for loop = 1:hyb


    path_to_fish = ['path=[' PathName '\' 'pos' num2str(posnum) '\' num2str(loop) '.tif' ']'];

    MIJ.run('Open...', path_to_fish);

    MIJ.run('Split Channels');
    for c = 1:length(channels)
        name = ['C' num2str(c) '-' num2str(loop) '.tif'];
        img(loop).color{c} = uint16(MIJ.getImage(name));
    end
    zleng = size(img(loop).color{1},3);
    MIJ.run('Close All');
    
    img(loop) = correctbackground(img(loop), corrections(loop),channels);
    
    for d = 1:length(channels)
        namesh{channels(d)} = ['C' num2str(channels(d)) '-'  num2str(loop) '.tif'];
        MIJ.createImage(namesh{channels(d)}, img(loop).color{channels(d)}, true);
    end
    
    str = [];
    for channum = 1:length(channels)
        str = [str ' image' num2str(channum) '=C' num2str(channum) '-' num2str(loop) '.tif'];
    end
    MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
    MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(length(channels)) ' slices=' num2str(zleng) ' frames=1 display=Grayscale']);
    MIJ.run('Subtract Background...', 'rolling=3 stack');
    MIJ.run('Split Channels');
    

    for dum = 1:length(channels)
        name = ['C' num2str(dum) '-Concatenated Stacks'];
        im = uint16(MIJ.getImage(name));
        hybnum(loop).color{dum} = im;
    end
    
    
    
    MIJ.run('Close All');

    
    for j = 1:length(channels)
            hybnum(loop).color{j} = imwarp(hybnum(loop).color{j},tforms{j},'OutputView',imref3d(size(hybnum(loop).color{2})));
    end

end
MIJ.exit


