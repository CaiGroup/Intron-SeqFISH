function [seeds, int, ints] = numseeds(PosList,dotlocations,hybs)

%load([PathName '\Pos' num2str(Posnum) '\pos' num2str(Posnum) 'Barcodes11092016.mat']);


seeds = PosList;
int = PosList;
ints = PosList;
for i = 1:size(dotlocations,2)
    for j = 1:size(dotlocations(i).cell,1)
        ind = find(strcmpi(dotlocations(i).cell{j,1},PosList));
        if size(PosList{ind,i+1},1) == 1
            seeds{ind,i+1} = size(dotlocations(i).cell{j,4},1);
            int{ind,i+1} = dotlocations(i).cell{j,7};
            ints{ind,i+1} = dotlocations(i).cell{j,6};
        else
            idx = dsearchn(PosList{ind,i+1},dotlocations(i).cell{j,2});
            holdint = zeros(size(PosList{ind,i+1},1),hybs);
            holdints = zeros(size(PosList{ind,i+1},1),hybs);
            for k = 1:size(PosList{ind,i+1},1)
                kee = idx ==k;
                hybnums = cell2mat(dotlocations(i).cell(j,4));
                hybnums = hybnums(kee);
                inttemps = cell2mat(dotlocations(i).cell(j,6));
                inttemp = cell2mat(dotlocations(i).cell(j,7));
                inttemps = inttemps(kee);
                inttemp = inttemp(kee);
                holdint(k,hybnums) = inttemp;
                holdints(k,hybnums) = inttemps;
            end
            holder = sum(holdints>0,2)';
            seeds{ind,i+1} = holder;
            int{ind,i+1} = holdint;
            ints{ind,i+1} = holdints;
        end
    end
end
%save([PathName '\Pos' num2str(Posnum) '\pos' num2str(Posnum) 'Barcodes11092016.mat'],'seeds','-append');