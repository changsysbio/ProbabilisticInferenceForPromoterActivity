function attribute = CombineAttribute(info, properties)

    % This function returns the combined attribute from lineageSpecific.
    %
    % Examples:
    % attribute = CombineAttribute(experimentInformation, {'cellSize'});
    % attribute = CombineAttribute(experimentInformation, ...
    %     {'cellSize' 'fluorescence'});

    attribute = [];
    for j = 1 : length(properties)
        eval(['attribute.' properties{j} ' = [];']);
    end

    for i = 1 : length(info.lineageSpecific)
        for j = 1 : length(properties)
            previous = getfield(attribute, properties{j});
            if size(previous, 1) > 1
                previous = previous';
            end

            current = getfield(info.lineageSpecific(i), properties{j});
            if size(current, 1) > 1
                current = current';
            end

            attribute = setfield(attribute, properties{j}, [previous current]);
        end
    end

end
