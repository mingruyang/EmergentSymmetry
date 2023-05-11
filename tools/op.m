function r = op(x)
    N = length(x);
    r = kron(x{1},x{2});
    for i = 3:N
        r = kron(r,x{i});
    end
end
