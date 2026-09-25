### Dataset description ###
# Representative c.m.-energy point of the 76-point scan (4.226–4.95 GeV, 15.67 fb^-1 total)
data_4226  = DatasetManager.real_data.find("703_4230")      # 4.226 GeV real data sample
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")   # inclusive MC at the same energy

# Signal decay card (EvtGen syntax): e+e- -> D_s*+ D_s*-, D_s* -> gamma D_s, D_s -> K+ K- pi
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s*-    PHSP;
    Enddecay

    Decay D_s*+
    1.0000 gamma D_s+     PHSP;
    Enddecay

    Decay D_s*-
    1.0000 gamma D_s-     PHSP;
    Enddecay

    Decay D_s+
    1.0000 K+ K- pi+      PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-      PHSP;
    Enddecay

    End
DECAYCARD

# 500k signal exclusive MC events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "dsstar_dsstar_signal_mc"
  config.related_dataset = data_4226
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :straight_line   # scan / line-shape aware cross section
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "DsStarDsStar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.226]})
            .set_alias({"std::vector<double>" => "Vdouble"})

# Semi-inclusive selection: reconstruct ONE D_s* per event; the other D_s* is the recoil.
event_selection = Selection.new
  .select_track {                       # charged-track quality cuts
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
  }
  .select_photon {                      # photon for the D_s* -> gamma D_s transition
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025             # E > 25 MeV (barrel)
    energyThreshold_e 0.025             # E > 25 MeV (endcap)
    nGam              ">=1"             # at least one radiative photon
  }
  .pid(method: :probability) {          # kaon / pion identification
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ / K-
    identify :pion, against: [:kaon, :proton]   # pi+ / pi-
    nkp  ">=1"                          # at least two kaons (K+ and K-)
    nkm  ">=1"
    npip ">=1"                          # at least one pion (pi+ or pi-)
    npim ">=1"
  }
  # Reconstruct one D_s* (recID 1 = D_s*+) from gamma D_s (D_s -> K+ K- pi) and
  # infer the opposite D_s* from the recoil four-momentum. This replaces the kinematic fit.
  .partial_rec([1]) do
    best_combination_by_mass :"D_s+", 1.9685   # build D_s candidate near nominal D_s mass
    require_recoil_mass 2.075, 2.150           # recoil (other D_s*) consistent with m_Ds*
  end

my_algorithm
  .note(:mass_window,
        "D_s candidates are built requiring |M(K+K-pi) - m_Ds| < 15 MeV; best_combination_by_mass keeps the closest combination and the explicit 15 MeV window is applied in the ROOT analysis.")
  .note(:missing_mass_formula,
        "Modified missing mass M_miss = m_miss + m(gamma K+K-pi) - m_Ds*; the m(gamma K+K-pi) - m_Ds* offset is not expressible in the DSL and is applied in the ROOT analysis.")
  .note(:efficiency_curve,
        "The missing-mass window |M_miss - m_Ds*| < 5 sigma_Mmiss(E_CM) is c.m.-energy dependent; sigma_Mmiss(E_CM) is parameterised per energy point, so only a fixed loose recoil window is applied in BOSS and tightened in ROOT.")
  .note(:signal_yield_fit,
        "The signal yield is extracted from a fit to the M(gamma K+K-pi) distribution (ROOT-level).")
  .note(:isr_vp_correction,
        "Born cross sections are corrected for ISR / vacuum-polarisation effects using an iterative line-shape treatment; this correction is applied outside BOSS.")
  .note(:energy_scan,
        "The measurement covers 76 c.m. energies from 4.226 to 4.95 GeV (15.67 fb^-1 total); the same algorithm and signal MC are applied at every energy point (4.226 GeV shown here as the representative sample).")
  .note(:charge_conjugate,
        "Both D_s*+ and D_s*- decays are generated; reconstructing one D_s* per event is applied identically to the charge-conjugate configurations.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([data_4226, incMC_4226, exMC_signal])