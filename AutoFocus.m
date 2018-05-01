function zstart = AutoFocus(hybnum, hyb, channels,debug) 

for loop = 1:hyb 
    zsize(loop) = size(hybnum(loop).color{1},3);
    for j = 1:length(channels)
        for i = 1:20 
            FM(i) = fmeasure(hybnum(loop).color{j}(:,:,i), 'LAPD',[128 128 256 256]);
        end
        if debug == 1
            figure;
            plot(FM)
        end
        dFM = diff(FM);
        [~,zpos(loop,j)] = max(FM);
    end
    zstart(loop) = round(mean(zpos(loop,:)));
end