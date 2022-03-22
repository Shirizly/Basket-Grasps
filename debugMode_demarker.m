if ~isempty(labels)
    for i = 1:numel(labels)
        delete(labels{i});
    end
end
if ~isempty(lines)
    for i = 1:numel(lines)
        delete(lines{i});
    end
end