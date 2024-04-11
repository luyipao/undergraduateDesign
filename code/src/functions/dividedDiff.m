function result = dividedDiff(x, y)
n = length(x);
diff_table = zeros(n, n);
diff_table(:, 1) = y(:);
for j = 2:n
    for i = j:n
        if x(i) == x(i-j+1)
            diff_table(i,j) = 0;
            continue;
        end
        diff_table(i, j) = (diff_table(i, j - 1) - diff_table(i - 1, j - 1)) / (x(i) - x(i - j + 1));
    end
end
result = diff_table(n, n);
end