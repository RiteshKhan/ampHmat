function varargout = chop_mat(varargin)

    global opt;

    if nargin == 1
        % ---------------------------------------------------------
        % Mode 1: M -> chop(M)
        % ---------------------------------------------------------
        M = varargin{1};
        varargout{1} = mchop(M);

    elseif nargin == 3
        % ---------------------------------------------------------
        % Mode 2: U,S,V -> chopped low-rank factors
        % ---------------------------------------------------------
        U = varargin{1};
        S = varargin{2};
        V = varargin{3};

        U = mchop(U);
        V = diag(S) * mchop(V');

        varargout{1} = U;
        varargout{2} = V;

    else
        error('Usage: chop_mat(M) or chop_mat(U,S,V).');
    end
end