function [coef,pos] = optranslator(pr, N)
    c = zeros(2,N);
    p = zeros(2,N);
    for i = 1:N
        if pr(i) == 'A'
            c(:,i) = [1/2; 1i/2];
            p(:,i) = [1; 2];
        elseif pr(i) == 'B'
            c(:,i) = [1/2; -1i/2];
            p(:,i) = [1; 2];
        elseif pr(i) == 'C'
            c(:,i) = [1/2; 1/2];
            p(:,i) = [0; 3];
        elseif pr(i) == 'D'
            c(:,i) = [1/2; -1/2];
            p(:,i) = [0; 3];
        end
    end
    coef = ones(2^N,1);
    pos = zeros(2^N,1);
    for j = 1:2^N
        ind = de2bi(j-1,N,'left-msb');
        for i = 1:N
            coef(j) = coef(j)*c(ind(i)+1,i);
            pos(j) = pos(j) + p(ind(i)+1,i)*4^(N-i);
        end
    end
end