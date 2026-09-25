# ============================================================================
# e+e- -> D*+ anti-D*0 pi-  at sqrt(s) = 4.260 GeV (charge conjugation implied)
# Partial-reconstruction analysis: bachelor pi- and D+ are tagged,
# the anti-D*0 is inferred from the D+ pi- recoil.
# ============================================================================

### Dataset description ###
data_4260  = DatasetManager.real_data.find("703_4260")      # 827 pb^-1 real data at 4.260 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")   # corresponding inclusive MC sample

# Decay card for the signal process e+e- -> D*+ anti-D*0 pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ anti-D*0 pi- PHSP;
    Enddecay

    Decay D*+
    1.000 D+ pi0 PHSP;
    Enddecay

    Decay anti-D*0
    1.000 anti-D0 pi0 PHSP;
    Enddecay

    Decay D+
    1.000 K- pi+ pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive PHSP MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_DstDstbarPi"
  config.related_dataset = data_4260
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "DstDstbarPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                        # charged-track selection
        cos_theta 0.93                     # |cos(theta)| < 0.93
        Vz        10.0                     # |Vz| < 10 cm
        Vr        1.0                      # Vr < 1 cm
        nChrp     ">=2"                    # at least 2 positive tracks
        nChrn     ">=2"                    # at least 2 negative tracks (>= 4 charged in total)
    }
    .select_photon {                       # photon selection
        tdc_emc_start     0                # TDC window 0-14
        tdc_emc_end       14
        energyThreshold_b 0.025            # E > 25 MeV in the barrel
        energyThreshold_e 0.050            # E > 50 MeV in the endcap
        angle_to_track    10.0             # at least 10 deg from any charged track
        nGam              ">=2"            # at least two photons
    }
    .pid(method: :probability) {           # probability PID with K/pi separation
        prob_cut 0.001                     # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]   # K+ / K-
        identify :pion, against: [:kaon, :proton]   # pi+ / pi-
        nkm  ">=1"                         # >= 1 K-
        npip ">=2"                         # >= 2 pi+
        npim ">=1"                         # >= 1 pi-
    }
    # D+ -> K- pi+ pi+ : secondary vertex fit, best vertex chi2
    .secondary_vertex_fit([:km, :pip, :pip]) {
        build_virtual_particle(:D_plus).by_minimizing_verfit_chi2
    }
    .invariant_mass_of(:km, :pip, :pip).between(1.854, 1.884)    # D+ mass window 1.854-1.884 GeV/c^2
    # soft pi0 -> gamma gamma : mass-constrained Kalman fit
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).within(0.120, 0.145)   # pi0 mass window 0.120-0.145 GeV/c^2
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                                              # chi2 < 25
        npi0 ">=1"                                               # at least one pi0
    }
    # Partial reconstruction: tag bachelor pi- (recID 3) and D+ (recID 4, daughters
    # auto-expanded); the anti-D*0 (recID 2) is inferred from the D+ pi- recoil.
    .partial_rec([3, 4]) {
        best_combination_by_mass :D_plus, 1.86962                # D+ candidate closest to nominal mass
        require_recoil_mass 2.135, 2.175                         # D+ pi- recoil mass window
    }

my_algorithm
    .note(:soft_pi0_selection, "soft pi0 (from D*+ -> D+ pi0) accepted if at least one pi0 satisfies either 2.008 < M(D+pi0)-M(D+)+m(D+)-M(pi0)+m(pi0) < 2.013 GeV/c^2 or P*(pi0) in (0.03, 0.05) GeV in the D+ pi- recoil system; not expressible in the DSL")
    .note(:background_veto, "two-body background suppression applied after the partial reconstruction: require RM(D+)+M(D+)-m(D+) > 2.3 GeV/c^2 and reject events with RM(pi-) > 4.1 GeV/c^2")
    .note(:combinatorial_background, "wrong-sign events (D+ pi+ instead of D+ pi-) scaled by a factor 1.9 are used to model the combinatorial background")
    .note(:zc_cascade, "resonant cascade e+e- -> Zc(4025)+ pi-, Zc(4025)+ -> D*+ anti-D*0 treated separately with its own decay card and exclusive MC")
    .note(:dst_background, "D** (excited D meson) backgrounds e+e- -> D** anti-D*0 pi- (and charge conjugates) treated separately with dedicated exclusive MC samples")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([data_4260, incMC_4260, exMC_signal])