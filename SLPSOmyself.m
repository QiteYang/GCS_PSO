function ArcDec = SLPSOmyself(ArcDec,PopDec,Xmean)
    [N,D]       = size(ArcDec);
    [Nbest,~]   = size(PopDec); 
    %% SLPSO
    for i = 1 : N
        Better = zeros(1,D);
        for j = 1 : D
            id = unidrnd(Nbest);
            Better(j) = PopDec(id,j);
        end
        r1 = rand(D,1); r1 = r1';
        r2 = rand(D,1); r2 = r2';
        for kk = 1 : D
            if rand(1) > 0.8
                r1(kk) = 0;
                r2(kk) = 0; 
            end
        end
        ArcDec(i,:) = ArcDec(i,:) + r1.*(Better-ArcDec(i,:)) + r2.*(Xmean-ArcDec(i,:));
    end
end