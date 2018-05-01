
PathName = uigetdir;

    
    listing = dir([PathName '\Pos*']);
    
    for i = 1:length(listing)
        num(i) = str2double(listing(i).name(4:end));
    end

    keep = num <30;

    num = num(keep);
    corrall = [];
%p = 11;
    for p = 1:length(num)
        posnum = num(p);
        disp('Processing Images....')
        [ hybnum] = PreprocessLinus( channels, hyb, PathName, posnum,tforms,corrections);

        disp('Registering Images....')
        %Register Images
        [hybnum, corr_offset, images] = xCorrRegLinus( PathName,hybnum,hyb,channels,posnum );

        save([PathName '\pos' num2str(posnum) '\hybnum' num2str(posnum) 'v1.mat'],'hybnum','corr_offset','images','-v7.3')
        %clear hybnum
        corrall = [corrall; corr_offset];
    end
    %save([PathName '\corrfocuspos11.mat'],'corrall','zfocusall')
    
 