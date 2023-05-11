function cs = Gtranslator(g, place, op, N)
    cs = zeros(2^(2*N),1);
    for n = 1:length(place)
        pr = op(place(n),:);
        [coef,pos] = optranslator(pr, N);
        for j = 1:2^N
            cs(pos(j)+1) = cs(pos(j)+1) + g(n)*coef(j);
        end
    end
end