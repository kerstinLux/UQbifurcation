function [numRows,M] = build_PowerMatrix(N,s)
% builds a matrix, where #columns=N (truncation of PCE) and
% #rows=#combinations of N summation indices

M=zeros(1,N);
% pure diagonal entries
for n=1:s-1
    M = [M; diag(n*ones(1,N))];
end

% diagonal entries with fixed columns of values 1,...,s-1-1 added
for j=1:s-1-1 % number of fixed columns to add
for n=1:s-1-j % value in fixed column to add
    C = combnk(1:N,j);
    for row = 1:length(C(:,1)) % row is number of possible combinations to distribute the j columns
        col = n*ones(N-1-C(row,end)+1,1);
        M_add = zeros(N-1-C(row,end)+1,N);
        M_add(:,C(row,:)) = repmat(col,1,length(C(row,:)));
        
        for i=1:s-1-n*j % i is a number that is allowed to add based on the already fixed value n, j times with upper bound s-1
            M_add_base = M_add;
            diag_add = diag(i*ones(1,N-C(row,end)));
            diagset = 1:N-C(row,end);
            if C(row,end)~=N
                M_add(:,C(row,end)+1:end) = diag_add(1:end,diagset);
            end
            M = [M; M_add];
            M_add = M_add_base;
        end
    end
end
end
numRows = length(M(:,1));
end