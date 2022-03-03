function z = zip(a,b)
% zip combines two equally sized arrays column wise such that
% [column 1 of a, column 1 of b, column 2 of a, column 1 of b, ...]

z = [zeros(size(a)),zeros(size(a))];
for idx = 1:size(a,2)
    z(:,2*idx-1) = a(:,idx);
    z(:,2*idx) = b(:,idx);
end
end