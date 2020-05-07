function F = interp2a(varargin)
% A modification of the interp2(a,b,c,d,e,'linear',0) to speed up

ExtrapVal = 0;
arg1= varargin{1};
arg2= varargin{2};
arg3= varargin{3};
arg4= varargin{4};
arg5= varargin{5};

% linear(extrapval,x,y,z,s,t), X and Y specified.
[nrows,ncols] = size(arg3);
mx = numel(arg1); my = numel(arg2);
if (mx ~= ncols || my ~= nrows) && ~isequal(size(arg1),size(arg2),size(arg3))
    error('MATLAB:interp2:linear:XYZLengthMismatch',...
        'The lengths of the X and Y vectors must match Z.');
end
if nrows < 2 || ncols < 2
    error('MATLAB:interp2:linear:sizeZ','Z must be at least 2-by-2.');
end
s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);

 
if nrows < 2 || ncols < 2
    error('MATLAB:interp2:linear:sizeZsq','Z must be at least 2-by-2.');
end
if ~isequal(size(s),size(t))
    error('MATLAB:interp2:linear:XIandYISizeMismatch',...
        'XI and YI must be the same size.');
end

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols));
if ~isempty(sout), s(sout) = 1; end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows));
if ~isempty(tout), t(tout) = 1; end

% Matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% Compute interpolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute interpolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Now interpolate.
onemt = 1-t;
F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
     ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;

% Now set out of range values to ExtrapVal.
if ~isempty(sout), F(sout) = ExtrapVal; end
if ~isempty(tout), F(tout) = ExtrapVal; end





