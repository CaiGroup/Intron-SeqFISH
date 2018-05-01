function [ foundbarcodes, points, copynumfinal, rawfound, hybnum, copynumfinalsum, totdropped ] = BarcodeNoMiji_v5( channels, hybnum, hyb, multiplier, HCRorFISH, barcodekey, PathName, posnum, segmentation, mask,mask2,radius,conthresh,alloweddiff,scaledonoff, debug )
 
if isempty(mask)
    mask = ones(size(hybnum(1).color{1},1),size(hybnum(1).color{1},2));
end
if isempty(mask2)
    mask2 = mask;
end

if isempty(PathName)
    PathName = uigetdir;
end
channum = length(channels);
disp('Finding All Dots....')
%find all dots in all images of all hybs
alloweddiff = alloweddiff + 1; 

for j = 1:hyb
    if size(hybnum(j).color{1},3) == size(mask,3)
        for i = 1:length(channels); 
            hybnum(j).color{i} = hybnum(j).color{i}.*uint16(mask); 
        end
    elseif size(hybnum(j).color{1},3) > size(mask,3)
        mask(:,:,size(mask,3):size(hybnum(j).color{1},3)) = repmat(mask(:,:,end),1,1,size(hybnum(j).color{1},3)-size(mask,3)+1);
        for i = 1:length(channels); 
            hybnum(j).color{i} = hybnum(j).color{i}.*uint16(mask); 
        end
    elseif size(hybnum(j).color{1},3) < size(mask,3)
        mask = mask(:,:,1:size(hybnum(j).color{1},3));
        for i = 1:length(channels); 
            hybnum(j).color{i} = hybnum(j).color{i}.*uint16(mask); 
        end
    end
end

fullpath = [PathName '\pos' num2str(posnum) '\RoiSet'];
vertex = selfseg(fullpath);
numcells = length(vertex);

for i = 1:hyb
    [points(i).m,points(i).dots] = findDotsBarcodeV2(hybnum(i).color, multiplier(i,:), HCRorFISH,debug); %250-350
end

%make mask2 same size
for j = 1:hyb
    if size(hybnum(j).color{1},3) == size(mask2,3)
        maskhyb(j).mask = mask2;
    elseif size(hybnum(j).color{1},3) > size(mask2,3)
        masktemp = mask2;
        masktemp(:,:,size(mask2,3):size(hybnum(j).color{1},3)) = repmat(mask2(:,:,end),1,1,size(hybnum(j).color{1},3)-size(mask2,3)+1);
        maskhyb(j).mask = masktemp;
    elseif size(hybnum(j).color{1},3) < size(mask2,3)
        masktemp = mask2;
        masktemp = mask2(:,:,1:size(hybnum(j).color{1},3));
        maskhyb(j).mask = masktemp;
    end
end

for i = 1:hyb
    for j = 1:length(channels)
        linIndrem = sub2ind(size(points(i).m{j}), points(i).dots(j).channels(:,2), points(i).dots(j).channels(:,1), points(i).dots(j).channels(:,3));
        tf = maskhyb(i).mask(linIndrem);
        points(i).dots(j).channels(~tf,:)= [];
        points(i).dots(j).intensity(~tf,:)= [];
        clear tf;
    end
end

for i = 1:hyb 
    for j = 1:length(channels)
        points(i).dots(j).scaledIntensity = double(points(i).dots(j).intensity)/mean(points(i).dots(j).intensity);
    end
end

disp('Initializing....')
%initialize barcode matrix
for i = 1:hyb
    for j = 1:channum
        c = size(points(i).dots(j).channels,1);
        filler = zeros(c,1);
        filler(:) = j; 
        idxfill = zeros(c,1);
        idxfill(:) = 1:c;
        foundbarcodes(i).found(j).channel(:,i) = num2cell(filler);
        foundbarcodes(i).found(j).idx(:,i) = num2cell(idxfill);
    end
end

disp('Colocalizing....')
%colocalize channels and output found codes and dots indices
for i = 1:hyb
    for k = 1:hyb
        if k ~= i
            for j = 1:channum
                calledall = [];
                calledidx = [];
                for l = 1:channum
                    [~, ~, called1, ~, pair ] = colocalizeBarcodeV2( hybnum(i).color{j},hybnum(k).color{l},points, i,k, [j l],radius,debug );
                    calledall = [calledall, l*double(called1)];
                    brat = zeros(length(called1),1);
                    if ~isempty(pair)
                        brat(pair(:,1)) = pair(:,2);
                    end
                    calledidx = [calledidx, brat];
                end
                calledall = mat2cell(calledall, ones(length(called1),1),channum);
                x = cellfun(@removezeros,calledall,'UniformOutput',0);
                %x = calledall;
                foundbarcodes(i).found(j).channel(:,k) = x;
                calledidx = mat2cell(calledidx, ones(length(called1),1),channum);
                %x = cellfun(@removezeros,calledidx,'UniformOutput',0);
                %x = calledidx;
                foundbarcodes(i).found(j).idx(:,k) = calledidx;
            end
        end
    end
end

rawfound = foundbarcodes;
disp('Determining Barcodes....')
%call Barcodes
totdropped = 0;
for k = 1:hyb
    for j = 1:channum
        drop = [];
        bobo = [];
        bratat = {};
        br = [];
        len = [];
        multi = [];
        rows = [];
     
        lentemp = size(foundbarcodes(k).found(j).channel,1);
        for i = 1:lentemp
            bobo = cell2mat(foundbarcodes(k).found(j).channel(i,:)); 
            bratat{i,1} = bobo; 
        end
        br = cellfun(@(x) sum(x == 0),bratat,'UniformOutput',0);
        
        %remove uncallable dots
        called = zeros(lentemp,1);
        drop=cellfun(@(x) x>1,br);
        foundbarcodes(k).found(j).channel(drop,:) ={0};
        foundbarcodes(k).found(j).idx(drop,:) ={0};
        bratat(drop) = {[0 0 0 0]};
        %reduce multiple matches
        br2 = cellfun(@removezeros2,bratat,'UniformOutput',0);
        len = cellfun(@length,br2);
        multi = len+cell2mat(br)>hyb;
        [rows] = find(multi ==1);
        for l = 1:length(rows)
            possbarcodes = [];
            Dposs = [];
            A = [];
            posscell = [];
            code = [];
            gene = [];
            possreal = [];
            vernoi = [];
            vernoiz = [];
            ind = [];
            set = [];
            vari = [];
            variz = [];
            possbarcodes = combvec(foundbarcodes(k).found(j).channel{rows(l),:})';
            Dposs = pdist2(possbarcodes,barcodekey.barcode,'hamming');
            A = Dposs == min(min(Dposs));
            [code, ~] = find(A == 1);
            if length(code) == 1 && Dposs(A)*hyb < alloweddiff
                posscell = num2cell(possbarcodes);
                foundbarcodes(k).found(j).channel(rows(l),:) = posscell(code,:);
                bratat{rows(l)} = cell2mat(posscell(code,:));
                %called(rows(l)) = gene;
            elseif Dposs(A)*hyb < alloweddiff & length(code) > 1
                possreal = possbarcodes(code,:);
                ind = foundbarcodes(k).found(j).idx(rows(l),:);
                caca = zeros(1,length(channels));
                caca(1,j) = ind{k};
                ind{k} = caca;
                for p = 1:size(possreal,1)
                    for h = 1:hyb
                        if possreal(p,h) > 0
                            set(p).set(h,:) = points(h).dots(possreal(p,h)).channels(ind{h}(possreal(p,h)),:);
                        else 
                            set(p).set(h,:) = zeros(1,3);
                        end
                    end
                    vari(p) = sum(var(set(p).set));
                    variz(p) = var(set(p).set(:,3));
                end
                vernoi = vari == min(vari);
                vernoiz = variz == min(variz);
                if sum(vernoi) == 1
                    foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(vernoi',:));
                    bratat{rows(l)} = possreal(vernoi',:);
%                 elseif sum(vernoi) > 1 && sum(vernoiz) == 1
%                     foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(vernoiz',:));
%                     bratat{rows(l)} = possreal(vernoiz',:);
                else
                    if scaledonoff == 1
                        for hybrid = 1:hyb
                            for numposs = 1:size(possreal,1)
                                if possreal(numposs,hybrid) == 0
                                    intensity(numposs,hybrid) = NaN;
                                else
                                    intensity(numposs,hybrid) = points(hybrid).dots(possreal(numposs,hybrid)).scaledIntensity(ind{hybrid}(possreal(numposs,hybrid)));
                                end
                            end
                        end
                        intdif = sum(abs(intensity-intensity(1,k)),2,'omitnan');
                        intensity = [];
                        if sum(intdif == min(intdif)) == 1
                            [~,I]=min(intdif);
                            foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(I,:));
                            bratat{rows(l)} = possreal(I,:);
                        else
                            totdropped = totdropped + 1;
                            foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                            foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                            bratat(rows) = {[0 0 0 0]};
                        end
                    else
                        totdropped = totdropped + 1;
                        foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                        foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                        bratat(rows) = {[0 0 0 0]};
                    end
                end
            else
                foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                bratat(rows) = {[0 0 0 0]};
%                 for n = 1:hyb
%                     mu = zeros(1,channum);
%                     mu(1,posscell{code,n}) = 1; 
%                     cell2{1,n} = mu;
%                 end
%                 foundbarcodes(k).found(j).idx(rows(l),:) = cellfun(@(x,y) x.*y,foundbarcodes(k).found(j).idx(rows(l),:),cell2,'UniformOutput',0);
            end
        end
        % Fill in dropped cells
        br2 = cellfun(@removezeros2,bratat,'UniformOutput',0);
        len = cellfun(@length,br2);
        missing = len == hyb-1;
        [rows] = find(missing ==1);
        for l = 1:length(rows)
            possbarcodes = [];
            Dposs = [];
            A = [];
            posscell = [];
            code = [];
            gene = [];
            possreal = [];
            vernoi = [];
            vernoiz = [];
            ind = [];
            set = [];
            vari = [];
            variz = [];
            possbarcodes = combvec(foundbarcodes(k).found(j).channel{rows(l),:})';
            Dposs = pdist2(possbarcodes,barcodekey.barcode,'hamming');
            A = Dposs == min(min(Dposs));
            [code, gene] = find(A == 1);
            if length(gene) == 1 && Dposs(A)*hyb < alloweddiff
                %posscell = num2cell(possbarcodes);
                foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(barcodekey.barcode(gene,:));
                bratat{rows(l)} = barcodekey.barcode(gene,:);
                %called(rows(l)) = gene;
           elseif Dposs(A)*hyb < alloweddiff & length(code) > 1
                possreal = possbarcodes(code,:);
                ind = foundbarcodes(k).found(j).idx(rows(l),:);
                caca = zeros(1,5);
                caca(1,j) = ind{k};
                ind{k} = caca;
                for p = 1:size(possreal,1)
                    for h = 1:hyb
                        if possreal(p,h) > 0
                            set(p).set(h,:) = points(h).dots(possreal(p,h)).channels(ind{h}(possreal(p,h)),:);
                        else 
                            set(p).set(h,:) = zeros(1,3);
                        end
                    end
                    vari(p) = sum(var(set(p).set));
                    variz(p) = var(set(p).set(:,3));
                end
                vernoi = vari == min(vari);
                vernoiz = variz == min(variz);
                if sum(vernoi) == 1
                    foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(vernoi',:));
                    bratat{rows(l)} = possreal(vernoi',:);
%                 elseif sum(vernoi) > 1 && sum(vernoiz) == 1
%                     foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(vernoiz',:));
%                     bratat{rows(l)} = possreal(vernoiz',:);
                else
                    if scaledonoff == 1
                        for hybrid = 1:hyb
                            for numposs = 1:size(possreal,1)
                                if possreal(numposs,hybrid) == 0
                                    intensity(numposs,hybrid) = NaN;
                                else
                                    intensity(numposs,hybrid) = points(hybrid).dots(possreal(numposs,hybrid)).scaledIntensity(ind{hybrid}(possreal(numposs,hybrid)));
                                end
                            end
                        end
                        intdif = sum(abs(intensity-intensity(1,k)),2,'omitnan');
                        intensity = [];
                        if sum(intdif == min(intdif)) == 1
                            [~,I]=min(intdif);
                            foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(I,:));
                            bratat{rows(l)} = possreal(I,:);
                        else
                            totdropped = totdropped + 1;
                            foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                            foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                            bratat(rows) = {[0 0 0 0]};
                        end
                    else
                        totdropped = totdropped + 1;
                        foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                        foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                        bratat(rows) = {[0 0 0 0]};
                    end
                end
            else
                foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                bratat(rows(l)) = {[0 0 0 0]};
%                 for n = 1:hyb
%                     mu = zeros(1,channum);
%                     mu(1,posscell{code,n}) = 1; 
%                     cell2{1,n} = mu;
%                 end
%                 foundbarcodes(k).found(j).idx(rows(l),:) = cellfun(@(x,y) x.*y,foundbarcodes(k).found(j).idx(rows(l),:),cell2,'UniformOutput',0);
            end
        end 
        % call barcodes
        if ~isempty(foundbarcodes(k).found(j).channel)
            cell2 = [];
            minmat = [];
            tester = [];
            logicalstuff = [];
            posscell = {};
            [posscell{1:size(foundbarcodes(k).found(j).channel,1),1:hyb}] = deal(0);
            %call dots
            D = pdist2(cell2mat(foundbarcodes(k).found(j).channel),barcodekey.barcode,'hamming');
            [r,c] = size(D);
            minmat = cellfun(@(x) x == min(x,[],2) & x<=((alloweddiff-1)/hyb),mat2cell(D,ones(1,r),c),'UniformOutput',0);
            %minmat = cellfun(@(x) x == 0,mat2cell(D,ones(1,r),c),'UniformOutput',0);
            [re, co] = cellfun(@(x) find(x),minmat,'UniformOutput',0);
            dasdf = cellfun(@isempty,re,'UniformOutput',0);
            re(logical(cell2mat(dasdf))) = {0};
            mm = cellfun(@length,re);
            mmdrop = mm > 1;
            re(mmdrop) = {0};
            rows = find(cell2mat(re));
            posscell(rows,:) = num2cell(barcodekey.barcode(cell2mat(co(rows)),:));
            co(logical(cell2mat(dasdf))) = {0};
            co(mmdrop) = {0};
            foundbarcodes(k).found(j).called = cell2mat(co);
            awe = foundbarcodes(k).found(j).called > 0;
            rawfound(k).found(j).called = foundbarcodes(k).found(j).called;
            for m = 1:size(foundbarcodes(k).found(j).channel,1)
                for n = 1:hyb
                    mu = zeros(1,channum);
                    if foundbarcodes(k).found(j).channel{m,n} > 0
                        mu(1,foundbarcodes(k).found(j).channel{m,n}) = 1;
                    end
                    cell2{m,n} = mu;
                end
            end
            foundbarcodes(k).found(j).idx = cellfun(@(x,y) x.*y,foundbarcodes(k).found(j).idx,cell2,'UniformOutput',0);
            %whereamI = sprintf('hyb %s, channel %s, pos %s',num2str(k),num2str(j),num2str(posnum));
            %disp(whereamI)
        else
            foundbarcodes(k).found(j).idx = [];
            %whereamI = sprintf('hyb %s, channel %s, pos %s',num2str(k),num2str(j),num2str(posnum));
            %disp(whereamI)
        end
    end
end
disp('Forming Point Consensus....')
%consensus point calling
for i = 1:hyb
    for j = 1:length(channels)
        foundbarcodes(i).found(j).compiled = zeros(size(foundbarcodes(i).found(j).called,1),hyb);
    end
end
% create a compiled list for every dot
disp('Compiling List....')
for i = 1:hyb
    for j = 1:length(channels)
        calledrows = find(foundbarcodes(i).found(j).called>0);
        if ~isempty(calledrows)
        [r,c,v]=cellfun(@(x) find(x),foundbarcodes(i).found(j).idx,'UniformOutput',0);
            for k = 1:length(calledrows)
                for l = 1:hyb
                    if ~isempty(c{calledrows(k),l})
                        foundbarcodes(l).found(c{calledrows(k),l}).compiled(v{calledrows(k),l},i) = foundbarcodes(i).found(j).called(calledrows(k));
                    end
                end
            end
        end
    end
end

disp('Dropping Ambiguous Matches....')
for i = 1:hyb
    for j = 1:length(channels)
        % Drop Ambigous Matches
        cellver = mat2cell(foundbarcodes(i).found(j).compiled,ones(size(foundbarcodes(i).found(j).compiled,1),1),hyb);
        cellver = cellfun(@removezeros, cellver,'UniformOutput',0);
        [M,F,C] = cellfun(@mode,cellver,'UniformOutput',0);
        eq = cellfun(@(x) length(x{1}),C);
        M = cell2mat(M);
        M(eq>1)=0;
        M(cell2mat(F)<conthresh) = 0;
        foundbarcodes(i).found(j).consensus = M;
    end
end

%only using hyb1 as seed
% for i = 1:hyb
%     for j = 1:length(channels)
%         foundbarcodes(i).found(j).consensus = foundbarcodes(1).found(j).compiled(:,1);
%     end
% end

disp('Assigning to Cell....')
if strcmp(segmentation, 'roi')
    %fullpath = [PathName '\pos' num2str(posnum) '\RoiSet'];
    %vertex = selfseg(fullpath);
    copynumfinal(:,:) = barcodekey.names;
    for i = 1:length(vertex)
        for j = 1:hyb
            %j = 1;
            allcalled = [];
            for k = 1:channum
                include(j).points(k).channel = inpolygon(points(j).dots(k).channels(:,1),points(j).dots(k).channels(:,2),vertex(i).x,vertex(i).y);
                allcalled = [allcalled; foundbarcodes(j).found(k).consensus(include(j).points(k).channel,:)];
            end
            copy = histc(allcalled(:),0:length(barcodekey.names));
            copynum(:,j) = copy(2:end);
        end
        copynumfinal(:,i+1) = num2cell(max(copynum,[],2));
        copynumfinalsum(:,i+1) = num2cell(sum(copynum,2));
    end
            
else
    copynumfinal(:,:) = barcodekey.names;

    for j = 1:hyb    
        allcalled = [];
        for k = 1:channum
            allcalled = [allcalled; foundbarcodes(j).found(k).consensus];
        end
        copy = histc(allcalled(:),0:length(barcodekey.names));
        copynum(:,j) = copy(2:end);
    end

    copynumfinal(:,2) = num2cell(max(copynum,[],2));
    copynumfinalsum(:,2) = num2cell(sum(copynum,2));
end

disp('Field Completed')