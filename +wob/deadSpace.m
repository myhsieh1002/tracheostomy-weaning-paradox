function ds = deadSpace(cfg, dev)
%DEADSPACE Total physiological dead space in place with a given device.
%
%   ds = wob.deadSpace(cfg, dev) returns a struct with the dead-space
%   budget in litres.
%
%   MODEL
%   -----
%   With a cuffed tube in situ, every breath enters through the tube lumen.
%   The airway above the tube TIP is therefore excluded from ventilation
%   and replaced by the tube's own internal volume:
%
%       V_D_total = (V_D_anatomic - V_D_bypassed) + V_D_apparatus
%
%   `V_D_bypassed` is the anatomic volume proximal to the tube tip. Both an
%   ETT and a tracheostomy tube are sited with the tip in the mid-trachea,
%   so both bypass essentially the SAME anatomic region -- mouth, pharynx,
%   larynx and upper trachea. What differs between them is the apparatus
%   volume replacing it: a ~27 cm ETT carries more internal volume than a
%   ~8 cm tracheostomy tube.
%
%   CONSEQUENCE, AND WHY IT MATTERS
%   -------------------------------
%   The often-quoted figure that "tracheostomy halves the dead space" is a
%   comparison against UNINSTRUMENTED breathing, not against an ETT. Once
%   the patient is already intubated, the upper-airway bypass term is
%   common to both arms and CANCELS in the ETT-vs-trach contrast. The
%   residual dead-space advantage of a tracheostomy is only the apparatus
%   difference, which is of order 5-10 mL against a tidal volume of ~450 mL.
%
%   The build spec anticipated a much larger dV_D as the mechanism behind
%   H3. Modelling it correctly makes H3's dead-space mechanism weaker than
%   assumed -- which does not undercut the programme thesis but sharpens
%   it: yet another axis on which the device moves less than expected.
%   `Vd_bypassed_mL` is exposed per device so this assumption is auditable
%   and can be varied in sensitivity analysis rather than buried.
%
%   Fields of `ds` (litres): anatomic, bypassed, apparatus, total,
%   anatomicEffective.
%
%   See also wob.simulateModeB, wob.requiredVentilation

arguments
    cfg (1,1) struct
    dev (1,1) struct
end

c = wob.constants();

anatomic_mL = cfg.patient.Vd_anatomic_mL;

% Per-device override, defaulting to the shared upper-airway value. Both
% device classes bypass the same region unless a config explicitly says
% otherwise.
if isfield(dev, 'Vd_bypassed_mL')
    bypassed_mL = dev.Vd_bypassed_mL;
else
    bypassed_mL = cfg.patient.Vd_upper_bypassed_mL;
end

if bypassed_mL > anatomic_mL
    error('wob:deadSpace:bypassExceedsAnatomic', ...
        ['Bypassed dead space (%.1f mL) exceeds total anatomic dead space (%.1f mL) ' ...
         'for device "%s". Check the config.'], bypassed_mL, anatomic_mL, dev.name);
end

apparatus_mL = dev.Vd_apparatus_mL;
anatomicEffective_mL = anatomic_mL - bypassed_mL;
total_mL = anatomicEffective_mL + apparatus_mL;

ds.anatomic          = anatomic_mL * c.ML_TO_L;
ds.bypassed          = bypassed_mL * c.ML_TO_L;
ds.apparatus         = apparatus_mL * c.ML_TO_L;
ds.anatomicEffective = anatomicEffective_mL * c.ML_TO_L;
ds.total             = total_mL * c.ML_TO_L;
end
