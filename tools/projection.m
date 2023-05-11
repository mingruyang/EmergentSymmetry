function O = projection(A, B)
    O = A - trace(B'*A)/trace(B'*B)*B;
end
