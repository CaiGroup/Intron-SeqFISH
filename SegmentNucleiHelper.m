

PathName = uigetdir;
listing = dir([PathName '\Pos*']);

for i = 0:length(listing)-1
%for i = 12:13
    fld = pwd;
    Miji;
    cd(fld);
    path_to_fish = ['path=[' PathName '\Pos' num2str(i) '\AllHybRegistrationCheck.tif' ']'];
    MIJ.run('Open...', path_to_fish);
    %MIJ.run('Z Project...', 'projection=[Max Intensity]');
    MIJ.run('Duplicate...', 'duplicate channels=1');
    MIJ.run('Save', ['save=[' PathName '\Pos' num2str(i) '\ClickDAPI.tif' ']']);
    MIJ.run('Close All');
    path_to_fish = ['path=[' PathName '\Pos' num2str(i) '\ClickDAPI.tif' ']'];
    MIJ.run('Open...', path_to_fish);
    MIJ.run('Z Project...', 'projection=[Max Intensity]');
    %MIJ.run('Enhance Contrast...', 'saturated=0.01 normalize equalize');
    MIJ.run('Gaussian Blur...', 'sigma=1');
    %MIJ.run('Auto Threshold', 'method=Mean');%didn't work well with 2i and
    %NIH3T3
    MIJ.run('Auto Threshold', 'method=Huang');%Initially used Li for 3T3 but some positions didn't work so switched to Huang
    MIJ.run('Convert to Mask');%activated
    %MIJ.run('Close')
    MIJ.run('Watershed')
    MIJ.run('Fill Holes');
    MIJ.run('Analyze Particles...', 'size=4000-Infinity pixel exclude clear add');
    MIJ.run('Close All')
    MIJ.run('Open...', path_to_fish);
    MIJ.run('Z Project...', 'projection=[Max Intensity]');
    pause('on');
    pause;
    MIJ.run('Close All')
    MIJ.exit;
end
