function [coef,place] = Gsignificant(G, N, tol)
    coef = [];
    place = [];
    for i = 1:2^N
        for j = 1:2^N
            if abs(G(i,j)) > tol
                coef(end+1) = G(i,j);
                place(end+1) = (j-1)*2^N+i;
            end
        end
    end
end