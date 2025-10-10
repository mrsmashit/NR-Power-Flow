function [A] = lowerfill(A)
    %% fills lower half of matrix
    A = triu(A).' + A - diag(diag(A));

end
