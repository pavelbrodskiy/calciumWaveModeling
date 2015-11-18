function matrix = processCellMatrix(cellMatrix)

matrixSize = size(cellMatrix);
fields = fieldnames(cellMatrix{1,1});

for i = 1:matrixSize(1)
    for j = 1:matrixSize(2)
        for k = 1:numel(fields)
            [i, j, k]
            matrix(i,j).(fields{k}) = median(double([cellMatrix{i, j}.(fields{k})]));
        end
    end
end

end