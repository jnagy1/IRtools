function k = subsasgn (k, index, val)
%  kronMatrix/subsasgn
%
%     changes the fields of a kronMatrix
%     
%     This only works for '.'
%

switch index.type
case '()'
    error('Array assignment not supported by kronmatrix objects')
case '{}'
    error('Cell array assignment not supported by kronmatrix objects')
case'.'
    if (index.subs == 'a')
        k.a = val;
    elseif (index.subs == 'b')
        k.b = val;
    else
        error ('Field Name unrecognized');
    end

end
