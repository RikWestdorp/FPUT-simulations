function [metadata, data] = load_simulation_output(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Failed to open file: %s', filename);
    end

    metadata = struct();
    data = struct();
    
    % Read metadata
    line = fgetl(fid);
    while ischar(line) && ~isempty(line)
        if startsWith(line, 'Metadata:')
            line = fgetl(fid);
            continue;
        end
        tokens = regexp(line, '(\w+)\s*=\s*([-\d\.eE+]+)', 'tokens');
        if ~isempty(tokens)
            key = tokens{1}{1};
            val = str2double(tokens{1}{2});
            metadata.(key) = val;
        end
        line = fgetl(fid);
    end

    % Read matrix blocks
    currentVar = '';
    matrix = [];
    while ischar(line)
        line = fgetl(fid);
        if ~ischar(line), break; end

        if isempty(line)
            continue;
        elseif endsWith(line, ':')
            % Save the previous variable if needed
            if ~isempty(currentVar)
                data.(currentVar) = matrix;
            end
            currentVar = strtrim(line(1:end-1));
            matrix = [];
        else
            nums = sscanf(line, '%f')';
            matrix = [matrix; nums];
        end
    end

    % Store the last variable
    if ~isempty(currentVar)
        data.(currentVar) = matrix;
    end

    fclose(fid);
end
