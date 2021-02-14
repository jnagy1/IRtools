classdef FreqMatrix
  %FREQMATRIX performs a frequency transformation
  %   Properties:
  %       transform: fft, dct, dwt
  %       imsize : vector containing the image size
  %       transpose: determines if the matrix is transposed
  %       wname: for dwt wavelet type (default is db1)
  %       wlevels: for dwt number of levels (default is 3)
  %
  % In matrix-vector multiplications,
  %       if input is vector, output is a vector
  %       if input is matrix, output is a matrix
  % J.Chung, updated 7/2017
  
  properties
    transform
    imsize
    transpose
    wname
    wlevels
  end
  
  methods
    
    function A = FreqMatrix (varargin)
      switch nargin
        case 0
          % default
          A.transform = '';
          A.imsize = '';
          A.transpose = 0;
          A.wname = '';
          A.wlevels = '';
        case 2
          A.transform = varargin{1};
          A.imsize = varargin{2};
          A.transpose = 0;
          A.wname = [];
          A.wlevels = '';
        case 3
          A.transform = varargin{1};
          A.imsize = varargin{2};
          A.wname = varargin{3};
          A.wlevels = '';
          A.transpose = 0;
        case 4
          A.transform = varargin{1};
          A.imsize = varargin{2};
          A.wname = varargin{3};
          A.wlevels = varargin{4};
          A.transpose = 0;
        otherwise
          error ('Incorrect number of arguments');
      end
      if isempty(A.wname)
        A.wname = 'db1';
      end
      if isempty(A.wlevels)
        A.wlevels = 3;
      end
    end
    
    function A = ctranspose(A)
      %  CTRANSPOSE The transpose of the FreqMatrix matrix
      if A.transpose == 0
        A.transpose = 1;
      else
        A.transpose = 0;
      end
    end
    
    function val = get(A, prop_name)
      %  Get FreqMatrix properties from specified object and return
      %  the value.
      switch prop_name
        case 'transform'
          val = A.transform;
        case 'imsize'
          val = A.imsize;
        case 'transpose'
          val = A.transpose;
        case 'wname'
          val = A.wname;
        case 'wlevels'
          val = A.wlevels;
        otherwise
          error([prop_name, 'Is not valid FreqMatrix property'])
      end
    end
    
    function A = set(A, varargin)
      % Set FreqMatrix properties and return the updated object.
      
      property_argin = varargin;
      
      while length( property_argin ) >= 2
        prop = property_argin{1};
        val = property_argin{2};
        property_argin = property_argin(3:end);
        
        switch prop
          case 'transform'
            A.transform = val;
          case 'imsize'
            A.imsize = val;
          case 'transpose'
            A.transpose = val;
          case 'wname'
            A.wname = val;
          case 'wlevels'
            A.wlevels = val;
          otherwise
            error('Valid FreqMatrix properties: transform, imsize, transpose, wname, wlevels')
        end
      end
      
    end
    
    function A = subsasgn(A, index, val)
      %  SUBSASGN  Define index assignment for transfromMatrix object.
      switch index.type
        case '()'
          error('Parenthetical indexing not supported for FreqMatrix object')
          
        case '.'
          switch index.subs
            case 'transform'
              P.transform = val;
            case 'imsize'
              P.imsize = val;
            case 'transpose'
              P.transpose = val;
            case 'wname'
              P.wname = val;
            case 'wlevels'
              P.wlevels = val;
            otherwise
              error('Invalid field names.');
          end
          
        case '{}'
          error('Cell array indexing not supported for FreqMatrix object.')
      end
      
    end
    
    
    function B = subsref(A, index)
      %  Define field name indexing for FreqMatrix object.
      switch index.type
        case '()'
          error('Paranthetical indexing not supported for FreqMatrix object.')
          
        case '.'
          switch index.subs
            case 'transform'
              B = A.transform;
            case 'imsize'
              B = A.imsize;
            case 'transpose'
              B = A.transpose;
            case 'wname'
              B = A.wname;
            case 'wlevels'
              B = A.wlevels;
            otherwise
              error('Invalid field name.')
          end
          
        case '{}'
          error('Cell array indexing not supported for FreqMatrix object.')
      end
    end
    
    function b = mtimes(A, x)
      %   Overload matrix multiplication operation for FreqMatrix
      % In matrix-vector multiplications,
      %       if input is vector, output is a vector
      %       if input is matrix, output is a matrix
      
      %   Note:  We scale fft so that the FreqMatrix is unitary.
      rs = 0;
      if isvector(x)
        x = reshape(x,A.imsize);
        rs = 1;
      end
      n = prod(A.imsize);
      
      if ~(isa(A, 'FreqMatrix'))
        error ('FreqMatrix Error: mtimes: Argument 1 must be a FreqMatrix')
      end
      
      switch A.transform
        case 'fft'
          if A.transpose
            b = sqrt(n) * ifftn(x);
          else
            b = fftn(x) / sqrt(n);
          end
        case 'dct'
          if A.transpose
            if isvector(x) %1D Inverse DCT
              b = idct(x);
            elseif ismatrix(x) %2D Inverse DCT
              b = idct2(x);
            end
          else
            if isvector(x) %1D DCT
              b = dct(x);
            elseif ismatrix(x) %2D DCT
              b = dct2(x);
            end
          end
        case 'dwt'
          if A.transpose
            if isvector(x) %1D Inverse Wavelet transform
              half = A.imsize(1)/2;
              cA = x(1:half);
              cD = x(half+1:end);
              b = idwt(cA, cD, A.wname);
            elseif ismatrix(x) %2D Inverse Wavelet transform
              [~,S] = wavedec2(x,A.wlevels,A.wname);
              b = waverec2(x(:),S,A.wname);
            end
          else
            if isvector(x) %1D Wavelet transform
              [cA,cD] = dwt(x,A.wname);
              b = [cA; cD];
            elseif ismatrix(x) %2D Inverse Wavelet transform
              [b,~] = wavedec2(x,A.wlevels,A.wname);
            end
          end
        otherwise
          error('Incorrect transform type.')
      end
      if rs
        b = b(:);
      else
        b = reshape(b,A.imsize);
      end
      
    end
    
    function varargout = size(A, dim) % Overload size
      d = [A.imsize(1)*A.imsize(2),A.imsize(1)*A.imsize(2)];
      
      if nargin == 2
        d = d(dim);
      end
      
      if nargout == 1 || nargout == 0
        varargout{1} = d;
      else
        for i = 1:length(d)
          varargout{i} = d(i);
        end
      end
      
    end
  end % methods
  
end

