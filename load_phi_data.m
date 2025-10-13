function [c_groups, c_values] = load_phi_data(filename)
    % Open and read all lines from the file
    fid = fopen(filename);
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = lines{1};

    % Remove header
    lines(1) = [];

    % Preallocate array
    num_rows = numel(lines);
    data = nan(num_rows, 9);

    % Process each line
    for i = 1:num_rows
        % Replace Fortran 'D' with 'E' for MATLAB compatibility
        clean_line = strrep(lines{i}, 'D', 'E');

        % Split into numeric tokens
        tokens = regexp(strtrim(clean_line), '\s+', 'split');

        % Convert to double
        data(i, :) = str2double(tokens);
    end

    % Determine block size and number of full blocks
    block_size = 1001;
    num_full_blocks = floor(num_rows / block_size);

    % Extract representative c values from each block
    c_column = data(:, 5);
    c_values = c_column(1:block_size:num_full_blocks * block_size);

    % Split data into blocks
    c_groups = mat2cell(data(1:num_full_blocks * block_size, :), ...
                        repmat(block_size, num_full_blocks, 1), 9);
end
