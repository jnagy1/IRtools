function Y = subsref(G, index)
%
%  Define field name indexing for kronMatrix object.
%
%               Y = subsref(A, index)
%
%          This is called whenever a subscripted reference to the
%          kronMatrix object is made.  Options are:
%              G.a
%              G.b
%              G.a{k}
%              G.b{k}
%          where k is an integer index.  That is, if A = G.a is a cell
%          array, then A{k} = G.a{k}.
%

%  J. Nagy  6/5/03

Y = G;

switch index(1).type

case '.'
  switch index(1).subs
  case 'a'
    Y = G.a;
  case 'b'
    Y = G.b;
  otherwise
    error('Subscript should be ''a'' or ''bb''')
  end

otherwise
  error('Array and cell array indexing are not supported.')
end

if length(index) == 2
  switch index(2).type
  case '{}'
    if isa(Y, 'cell')
      Y = Y{index(2).subs{1}};
    else
      if index(2).subs{1} == 1
        return
      else
        error('Index exceeds cell array dimension.')
      end
    end
  otherwise
    error('Only cell indexing is supported for factors.')
  end
elseif length(index) > 2
  error('Illegal reference.')
end

if length(Y) == 1
  Y = Y{1};
end
 




