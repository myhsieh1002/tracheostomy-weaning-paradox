function R = deviceResistance(dev, flow)
%DEVICERESISTANCE Flow-dependent resistance of an airway device.
%
%   R = wob.deviceResistance(dev, flow) returns K1 + K2*|flow| in
%   cmH2O/(L/s). This is the secant resistance dP/flow, not the tangent
%   (differential) resistance dP/dflow, which would be K1 + 2*K2*|flow|.
%
%   See also wob.rohrerDrop

arguments
    dev  (1,1) struct
    flow double
end

R = dev.K1 + dev.K2 .* abs(flow);
end
