function [ outputFlag, fileName, p ] = arginParser( p, size, argumentsPassed )
% Note to future-Pavel PLEASE PLEASE switch to input parser instead of eval

outputFlag = 0;

if mod(length(argumentsPassed),2) ~= 0
    disp('Inputs should be in the form (''Variable name'', Parameter Value)');
    outputFlag = 4;
    return
end

for i = 1:2:length(argumentsPassed)
    temp = argumentsPassed{i};
    if ~(temp(1)=='$') % If the input is not compartment-specific
        eval(['p.' temp ' = ' num2str(argumentsPassed{i+1}) ';']);
    else % If input is compartment-specific, this is indicated with $
        compartment = ones(size(1), size(2));
        
        if (temp(2)=='P') % Check if P or A
            compartment(:, 1:round(size(2)/2)) = 0;
        elseif (temp(2)=='A')
            compartment(:, round(size(2)/2):end) = 0;
        else
            error('$ flag used without A or P designation');
        end
        
        % Implement the change only in that compartment
        par = ['p.' temp(3:end)];
        if (isfield(p,temp(3:end)))
            eval([par ' = compartment .* ' num2str(argumentsPassed{i+1}) ' + ~compartment .* ' par ';']);
        else
            disp(['Parameter does not exist: ' par]);
            outputFlag = 5;
            return
        end
    end
end

fileName = '';
for i = 1:2:length(argumentsPassed)
    if ~(strcmp(argumentsPassed{i},'outputDirectory')||strcmp(argumentsPassed{i},'outputStart'))
        fileName = [fileName argumentsPassed{i} ' = ' num2str(argumentsPassed{i+1},'%.6f') ' '];
    end
end

end

