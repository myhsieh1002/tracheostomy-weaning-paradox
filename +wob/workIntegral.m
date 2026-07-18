function W = workIntegral(t, pressure, flow, mask)
%WORKINTEGRAL Mechanical work as the integral of pressure against volume.
%
%   W = wob.workIntegral(t, pressure, flow, mask) evaluates
%
%       W = int P dV = int P * flow dt
%
%   over the samples selected by `mask`, returning work in cmH2O*L.
%   Multiply by wob.constants().CMH2O_L_TO_JOULE for joules.
%
%   Expressing the work as int P*flow dt rather than int P dV lets a single
%   quadrature rule serve every term of the decomposition, so the parts are
%   guaranteed to sum to the total.
%
%   See also wob.ptpIntegral, wob.simulateModeB

arguments
    t        (:,1) double
    pressure (:,1) double
    flow     (:,1) double
    mask     (:,1) logical
end

W = trapz(t(mask), pressure(mask) .* flow(mask));
end
