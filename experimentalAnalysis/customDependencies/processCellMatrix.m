function matrix = processCellMatrix(cellMatrix)
% Processes the cell matrix of structures to obtain medians of the bins.

if (sum(sum(~cellfun(@isempty,cellMatrix))))
    tempMat = cellMatrix(~cellfun(@isempty,cellMatrix));
    matrixSize = size(cellMatrix);
    fields = fieldnames(tempMat{1,1});

    for i = 1:matrixSize(1)
        for j = 1:matrixSize(2)
            tempStruct = tempMat{1};
    matrix{i,j} = tempStruct(1);
        end
    end
    
for i = 1:matrixSize(1)
    for j = 1:matrixSize(2)
        for k = 1:numel(fields)
            if ~isempty(cellMatrix{i, j})
                matrix{i,j}.(fields{k}) = median(double([cellMatrix{i, j}.(fields{k})]));
            else
                matrix(i,j) = {[]};
            end
        end
    end
end

else
    
    matrix = [];
    
end

end