function [forest,feature] = ForestTrain(PopDec,PopObj,K,Type,FeaN)
[N,D] = size(PopDec);
forest = cell(1,K);
feature = cell(1,K);
nTrain = 4*ceil(N/5);
for i = 1 : K
    randIndex = randperm(N);
    randD = randperm(D);
    feature{i} = randD(1:ceil(FeaN*D));
    trainX = PopDec(randIndex(1:nTrain),feature{i});
    trainY = PopObj(randIndex(1:nTrain));
    if Type == 1   % classification tree for GLoc
        tree = fitctree(trainX,trainY);
    else     % regression tree for GCPD
        tree = fitrtree(trainX,trainY);
    end
    forest{i} = tree;
end
end