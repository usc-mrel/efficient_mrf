function parsave(fname, varargin)
    for ii = 1:2:length(varargin)
        val = varargin{ii};
        eval([varargin{ii+1} '= val;']);

        if ii == 1
            save(fname, varargin{ii + 1}, '-v7.3');
        else
            save(fname, varargin{ii + 1}, '-append');
        end
    end