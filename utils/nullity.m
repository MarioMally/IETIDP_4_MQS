function null = nullity(mat)
    [R] = qr(mat);
    [rows,~,~] = find(R);
    null = min(size(mat))-max(rows);
end