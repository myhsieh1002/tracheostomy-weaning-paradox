function PTP = ptpIntegral(t, pressure, mask)
%PTPINTEGRAL Pressure-time product over the selected samples.
%
%   PTP = wob.ptpIntegral(t, pressure, mask) evaluates
%
%       PTP = int P dt
%
%   returning cmH2O*s. Unlike work, the PTP does not vanish when flow is
%   zero, so it still counts isometric effort -- which is why it tracks the
%   oxygen cost of breathing better than WOB does.
%
%   See also wob.workIntegral

arguments
    t        (:,1) double
    pressure (:,1) double
    mask     (:,1) logical
end

PTP = trapz(t(mask), pressure(mask));
end
