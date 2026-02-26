function [XX,YY] = rotation(XX, YY, alpha)
% Function to rotate the curve of an angle alpha (centre of ration (0,0))

mat_rot=[cos(alpha)   -sin(alpha);
        sin(alpha)    cos(alpha)];

xy_XXYY_rot=mat_rot*[XX; YY];
XX=xy_XXYY_rot(1,:);
YY=xy_XXYY_rot(2,:);

end

