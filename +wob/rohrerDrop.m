function dP = rohrerDrop(dev, flow)
%ROHRERDROP Pressure drop across an airway device (Rohrer form).
%
%   dP = wob.rohrerDrop(dev, flow) returns the pressure drop in cmH2O for
%   flow in L/s.
%
%       dP = K1*flow + K2*flow*|flow|
%
%   The |flow| in the quadratic term keeps the drop odd-symmetric, so the
%   sign is correct on expiration as well as inspiration. The equivalent
%   flow-dependent resistance is dP/flow = K1 + K2*|flow|.
%
%   See also wob.deviceResistance, wob.calibrateRohrer

arguments
    dev  (1,1) struct
    flow double
end

dP = dev.K1 .* flow + dev.K2 .* flow .* abs(flow);
end
