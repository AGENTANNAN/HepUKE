# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# Real data
data_3650  = DatasetManager.real_data.find("709_3650")   # 3.650 GeV real data (44.5 pb^-1)
data_3773  = DatasetManager.real_data.find("712_3773")   # 3.773 GeV real data (2917 pb^-1)

# psi(3770) line-shape scan points covering 3.74-3.90 GeV
scan_data = %w[
  712_3768 712_3780 712_3800 712_3805 712_3815 712_3820
  712_3825 712_3830 712_3835 712_3840 712_3845 712_3855
  712_3860 712_3875 712_3880 712_3885
].map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")  # inclusive MC @ 3.650 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # inclusive MC @ 3.773 GeV

# Signal decay card: e+e- -> p pbar (no intermediate resonance; top mother psi(4260)
# is the BESIII/KKMC convention for continuum production)
decay_card_ppbar = <<~DECAYCARD
    Decay psi(4260)
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# ISR background decay card: psi(4260) -> gamma psi(3686), psi(3686) -> p pbar
decay_card_isr = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma psi(3686) PHSP;
    Enddecay

    Decay psi(3686)
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC for e+e- -> p pbar at every centre-of-mass energy point
# (3.650, 3.773 and every scan point) -- same card, one MC per energy point.
signal_points = [data_3650, data_3773] + scan_data
exMC_signal = DatasetManager.create_exclusive_mc_for(signal_points) do |config|
  config.sample_name   = "exmc_ppbar"     # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

# 500k-event ISR background sample associated with the 3.773 GeV data
exMC_isr = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_isr_psip_ppbar"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_isr
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "ppbar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})  # re-run per scan point with the corresponding CM energy
            .note(:background_veto,
                  "additional back-to-back cut requiring the CM proton-antiproton opening " \
                  "angle > 179 deg; no dedicated DSL implementation, applied as an extra cut " \
                  "outside the kinematic fit")

event_selection = Selection.new
event_selection.select_track {
        cos_theta 0.8      # |cos(theta)| < 0.8
        Vz        10.0     # |Vz| < 10 cm
        Vr        1.0      # Vr < 1 cm in the transverse plane
        nChrp     "==1"    # exactly one positively charged track
        nChrn     "==1"    # exactly one negatively charged track
        nNet      "==0"    # net charge zero
    }
    # No photon selection is applied.
    .pid(method: :probability) {
        prob_cut 0.001                              # PID probability > 0.001
        identify :proton, against: [:kaon, :pion]   # separate p+ / anti-p- from K and pi (both charges)
        nprp "==1"                                  # exactly one proton
        nprm "==1"                                  # exactly one antiproton
    }
    # Momentum window: |p - p_exp| < 40 MeV/c, p_exp = sqrt(s/4 - m_p^2)
    # (= 1.637 GeV/c at sqrt(s) = 3.773 GeV); tracks outside the window are dropped.
    .remove(:prp) { condition "abs(three_momentum_of(:prp) - sqrt(ECMS*ECMS/4.0 - 0.938272*0.938272)) > 0.040" }
    .remove(:prm) { condition "abs(three_momentum_of(:prm) - sqrt(ECMS*ECMS/4.0 - 0.938272*0.938272)) > 0.040" }
    # 4C kinematic fit constraining the p pbar system to the e+e- CM energy
    .kinematic_fit([:prp, :prm]) {
        nominal                  # nominal fit: corrected four-momenta are stored
        constrain_four_momentum  # 4C energy-momentum constraint
        chi2_cut 200             # loose chi2 cut; optimal cut applied later in ROOT
    }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_ppbar).apply(event_selection)

# Execute on real data, inclusive MC and exclusive MC samples
root_files = my_algorithm.execute_on(
  [data_3650, data_3773] + scan_data +
  [incMC_3650, incMC_3773] +
  exMC_signal + [exMC_isr]
)