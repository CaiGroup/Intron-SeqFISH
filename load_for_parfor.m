function hybnum =  load_for_parfor(PathName, i)

S = load([PathName '\' 'pos' num2str(i) '\' 'hybnum' num2str(i) 'v1.mat'],'hybnum');
hybnum = S.hybnum;