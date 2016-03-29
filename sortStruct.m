%% Function for sorting structs
% Ref: http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/
function out = sortStruct(features)
    headers = fieldnames(features);
    Acell = struct2cell(features);
    sz = size(Acell);
    % Convert to a matrix
    Acell = reshape(Acell, sz(1), []);      % Px(MxN)
    % Make each field a column
    Acell = Acell';                         % (MxN)xP
    % Sort by third field "v"
    Acell = sortrows(Acell, 3);
    % Put back into original cell array format
    Acell = reshape(Acell', sz);
    % Convert to Struct
    out = cell2struct(Acell, headers, 1);
    out = fliplr(out);
end