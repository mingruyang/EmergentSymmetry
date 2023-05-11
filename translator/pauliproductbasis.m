function op = pauliproductbasis(N)
    op = char(zeros(2^(2*N),N));
    for i = 1:2^N
        for j = 1:2^N
            ij = (j-1)*2^N+i;
            row = de2bi(i-1,N,'left-msb');
            col = de2bi(j-1,N,'left-msb');
            for k = 1:N
                if row(k) == 0 && col(k) == 1
                    op(ij,k) = 'A';
                elseif row(k) == 1 && col(k) == 0
                    op(ij,k) = 'B';
                elseif row(k) == 0 && col(k) == 0
                    op(ij,k) = 'C';
                elseif row(k) == 1 && col(k) == 1
                    op(ij,k) = 'D';
                end
            end
        end
    end
end
