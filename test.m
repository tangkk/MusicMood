root = 'testcase/4/';
subRoots = dir(root);
subRoots = subRoots(3:end); %exclude . and .. folder
numSongs = length(subRoots);
names = [];
datas = [];
for i = 1:1:numSongs
    name = subRoots(i).name;
    path = [root name];
    display(path);
    data = musicmood(path);
    names = [names name];
    datas = [datas data];
end

datas = reshape(datas, [6, length(datas)/6]);
datas = datas';