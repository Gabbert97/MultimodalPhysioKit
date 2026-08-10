function cell_matrix = string_to_cell_matrix_conversion(string_matrix)
%STRING_TO_CELL_MATRIX_CONVERSION Convert a string matrix to mixed cell data.
%
% The first column is preserved as text. Remaining columns are converted to
% numeric values with STR2DOUBLE.

cell_matrix = {};
[r,c] = size(string_matrix);
cell_matrix = cellstr(string_matrix(:,1));
for i = 1:r
    for j = 1:c-1
        cell_matrix{i,j+1} = str2double(string_matrix(i,j+1));
    end
end
end
