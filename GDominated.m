function [isDom,DomNo] = GDominated(offLoc,GPF)
isDom = false;
[N,~] = size(GPF);
DomNo = false(1,N);
for i = 1 : N
    k = any(offLoc<GPF(i,:)) - any(offLoc>GPF(i,:));
    if k > 0
        isDom = true;
    end
    if k > 0
        DomNo(i) = true;
    end
end
end