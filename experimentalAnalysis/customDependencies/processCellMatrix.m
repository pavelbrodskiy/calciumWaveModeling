function matrix = processCellMatrix(cellMatrix)
% Processes the cell matrix of structures to obtain medians of the bins.

matrixSize = size(cellMatrix);
fields = fieldnames(cellMatrix{1,1});

for i = 1:matrixSize(1)
    for j = 1:matrixSize(2)
        for k = 1:numel(fields)
            if cellMatrix{i, j}.flag
                matrix(i,j).(fields{k}) = median(double([cellMatrix{i, j}.(fields{k})]));
            else
                matrix(i,j).flag = median(double([cellMatrix{i, j}.flag]));
            end
        end
    end
end

end