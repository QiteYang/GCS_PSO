function y = ForestPredict(forest,decs,K,type,feature)
[N,~] = size(decs);
predY = inf(N,K);
for i = 1 : K
    tree = forest{i};
    feaD = feature{i};
    predY(:,i) = predict(tree,decs(:,feaD));
end
if type == 1 % classification for GLoc
    y = mode(predY,2);
else      % regression for GCPD
    y = mean(predY,2);
end
end