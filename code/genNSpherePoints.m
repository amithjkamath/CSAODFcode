function points = genNSpherePoints(n, varargin)
    if nargin == 1
        r = 1;
    else
        assert(nargin <= 2)
        r = varargin{1};
    end
    assert(n > 0)
    theta = 2.*pi.*rand([n 1]);
    phi   = acos(1 - 2.*rand([n 1]));
    x     = r.*sin(phi).*cos(theta);
    y     = r.*sin(phi).*sin(theta);
    z     = r.*cos(phi);
    points = [x y z];
end