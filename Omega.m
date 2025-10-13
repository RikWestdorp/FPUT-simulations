function Omega = Omega(f1, f2, g1, g2)
    % Assumes f1, f2, g1, g2 are all [T x N]
    
    % Cumulative sum along the second dimension (columns)
    cumsum_f2 = cumsum(f2, 2);  % [T x N]
    
    % Shift f1 one column to the right, insert 0 in the first column
    f1_shifted = [zeros(size(f1, 1), 1), f1(:, 1:end-1)];  % [T x N]
    cumsum_f1_shifted = cumsum(f1_shifted, 2);            % [T x N]
    
    % Element-wise multiplication and sum across columns (dimension 2)
    term1 = sum(cumsum_f2 .* g1, 2);       % [T x 1]
    term2 = sum(cumsum_f1_shifted .* g2, 2);  % [T x 1]
    
    % Final result: [T x 1]
    Omega = term1 + term2;  % still [T x 1]
end
