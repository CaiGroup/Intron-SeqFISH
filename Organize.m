function Organize( hyb, range )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
folder = uigetdir;

for i = hyb(1):hyb(2)
    listing(i).listing = dir([folder '/' num2str(i) '/*.tif']);
end



for i = 0:length(listing(hyb(1)).listing)-1
    [a, b] = regexp(listing(hyb(1)).listing(i+1).name,'\d*');
    roi = ['Pos' listing(hyb(1)).listing(i+1).name(a(end):b(end))];
    for j = hyb(1):hyb(2)
        file = strfind({listing(j).listing.name}, roi);
        pick = ~cellfun(@isempty,file);
        filename = listing(j).listing(pick).name;
        source = [folder '/' num2str(j) '/' filename];
        dest = [folder '/' range '/' roi '/' num2str(j) '.tif'];
        copyfile(source, dest);
    end
end

end

