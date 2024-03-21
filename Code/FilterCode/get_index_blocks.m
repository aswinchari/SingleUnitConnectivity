function blocks = get_index_blocks(v)
% blocks = get_index_blocks(v)
% Find blocks of consecutive integers in v and returns them in a cell
% array.

% Transpose input to a column vector
if size(v,1) == 1
    v = v';
end

% Sort input vector
v = sort(v);

% Find jumps
big_diff = diff(v)~=1;
big_diff = [0;big_diff];
jumps = find(big_diff);

% Fill blocks
jumps = [1;jumps];

for i = 1:length(jumps)
    if i < length(jumps)
        blocks{i} = v(jumps(i):(jumps(i+1)-1));
    else
        blocks{i} = v(jumps(i):end);
    end
end