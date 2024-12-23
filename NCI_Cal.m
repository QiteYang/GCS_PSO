function NCI = NCI_Cal(GLoc,i,Problem)
GD = sum(abs(GLoc - repmat(GLoc(i,:),Problem.N,1)),2);
neighbor = find(GD<=Problem.M);
neighborGLoc = GLoc(neighbor,:);
neighborGD = GD(neighbor);
neighborN = length(neighbor);
K = [];
for j = 1 : neighborN
    k = any(neighborGLoc(j,:)<GLoc(i,:)) - any(neighborGLoc(j,:)>GLoc(i,:));
    K = [K,k];
end
neighbor1 = find(K < 0);
neighbor2 = find(K == 0);
neighbor3 = find(K > 0);
neighbor1GD = neighborGD(neighbor1);
neighbor2GD = neighborGD(neighbor2);
neighbor3GD = neighborGD(neighbor3);
NCI = (length(neighbor1)+1)/(length(neighbor3)+1)*sum(Problem.M - neighbor1GD) + sum(Problem.M - neighbor2GD) + (length(neighbor3)+1)/(length(neighbor1)+1)*sum(Problem.M - neighbor3GD);
end