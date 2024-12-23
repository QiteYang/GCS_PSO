function choose = GridSelection(GLoc,K)
    [N,M] = size(GLoc);
    %% Calculate GR, GCD, and GD values of each solution
    GR   = sum(GLoc,2);
    GCD  = zeros(1,N);
    GD   = inf(N);
    for i = 1 : N-1
        for j = i+1 : N
            GD(i,j) = sum(abs(GLoc(i,:)-GLoc(j,:)));
            GD(j,i) = GD(i,j);
        end
    end
    %% Detect the grid-based dominance relation of each two solutions
    G = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(GLoc(i,:)<GLoc(j,:))-any(GLoc(i,:)>GLoc(j,:));
            if k == 1
                G(i,j) = true;
            elseif k == -1
                G(j,i) = true;
            end
        end
    end
    %% Environmental selection
    Remain = true(1,N);
    while sum(Remain) > N-K      
        % Choose the best one among the remaining solutions in the front
        CanBeChoose = find(Remain);
        temp  = find(GR(CanBeChoose)==min(GR(CanBeChoose)));
        temp2 = find(GCD(CanBeChoose(temp))==min(GCD(CanBeChoose(temp))));
        q = unidrnd(length(temp2));
        %[~,q] = min(Front(CanBeChoose(temp(temp2))));
        q     = CanBeChoose(temp(temp2(q)));
        Remain(q)   = false;      
        % Update the GCD values
        GCD = GCD+max(M-GD(q,:),0);
        % Update the GR values
        Eq  = GD(q,:)==0 & Remain;
        Gq  = G(q,:) & Remain;
        NGq = Remain.*(1-Gq);
        Nq  = GD(q,:)<M & Remain;
        GR(Eq) = GR(Eq) + M+2;
        GR(Gq) = GR(Gq) + M;
        PD = zeros(N,1);
        for p = find((Nq.*NGq).*(1-Eq))
            if PD(p) < M-GD(q,p)
                PD(p) = M - GD(q,p);
                Gp = G(p,:) & Remain;
                for r = find(Gp.*(1-(Gq+Eq)))
                    if PD(r) < PD(p)
                        PD(r) = PD(p);
                    end
                end
            end
        end
        pp = logical(NGq.*(1-Eq));       
        GR(pp) = GR(pp) + PD(pp);
    end
    choose = ~Remain;
end