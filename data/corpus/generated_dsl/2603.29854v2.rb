# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# Nine ψ(2S) energy-scan points (BOSS version 704), 3.5815–3.7098 GeV, 495 pb⁻¹ total.
# Sample name convention: [BOSS version]_[CMS energy in MeV]
scan_points = [
  DatasetManager.real_data.find("704_3581"),   # 3581.5 MeV
  DatasetManager.real_data.find("704_3670"),   # 3670.2 MeV
  DatasetManager.real_data.find("704_3680"),   # 3680.1 MeV
  DatasetManager.real_data.find("704_3682"),   # 3682.8 MeV
  DatasetManager.real_data.find("704_3684"),   # 3684.2 MeV
  DatasetManager.real_data.find("704_3685"),   # 3685.3 MeV
  DatasetManager.real_data.find("704_3686"),   # 3686.5 MeV
  DatasetManager.real_data.find("704_3691"),   # 3691.4 MeV
  DatasetManager.real_data.find("704_3709")    # 3709.8 MeV
]

# ConExc decay card for the continuum process e+e- -> K+K- (ConExc mode 45, K+K- PHSP).
# ConExc models ISR (up to 2nd order) and vacuum polarisation; the DSL auto-detects the
# literal token "ConExc", switches to the no-KKMC template and injects
# "Particle vpho <ECMS> 0.0" per scan energy point -- so "Particle vpho" is omitted here.
decay_card_kk = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 45;
    Enddecay
    End
DECAYCARD

# 100k-event exclusive signal MC at each of the nine energy points (same card / count,
# differing only by the related real dataset).
exmc_kk = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_kk_conexc"   # per-point suffix is appended automatically
  config.events        = 100_000
  config.decay_card    = decay_card_kk
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "ScanKpKm"
my_alg = Algorithm.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})   # per-point vpho energy handled by ConExc injection

event_selection = Selection.new
event_selection
  .select_track {          # Charged track selection
    cos_theta 0.93         # |cos(theta)| < 0.93
    Vz        10.0         # |Vz| < 10 cm
    Vr        1.0          # Vr < 1 cm in transverse plane
    nChrp     "==1"        # exactly one positive track
    nChrn     "==1"        # exactly one negative track
    nNet      "==0"        # net charge zero
  }
  .pid(method: :probability) {   # PID by the probability method
    prob_cut 0.001               # PID probability > 0.001
    identify :kaon, against: [:pion]   # K+ and K- (charge-conjugation shorthand), against pions
    nkp "==1"                    # one K+
    nkm "==1"                    # one K-
  }
  # 4C kinematic fit to K+K-, nominal, loose chi2 cut (tight cut applied in ROOT)
  .kinematic_fit([:kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

# BOSS-side procedures that the DSL cannot express, preserved as notes
my_alg
  .note(:bhabha_veto, "Bhabha (e+e- -> e+e-) background suppressed by requiring E/p < 0.8 for both
    charged tracks; E/p from EMC cluster energy over MDC momentum")
  .note(:cosmic_ray_veto, "cosmic-ray events suppressed by requiring |T(K+) - T(K-)| < 3 ns, where
    T is the TOF-measured flight time of each charged track")
  .note(:extra_track_veto, "events with extra charged tracks suppressed by requiring
    |p_CMS - p(K+) - p(K-)| < 0.07 GeV/c, where p_CMS is the expected total three-momentum
    of the K+K- system in the c.m. frame")
  .with_decay_card(decay_card_kk)
  .apply(event_selection)

# Run on the nine real-data scan points together with their exclusive signal MC.
# Signal-region window (M(K+K-) within 0.12 GeV of sqrt(s)) and the ISR / beam-energy-spread
# convoluted line-shape fit are applied in the ROOT-level unbinned maximum-likelihood fit.
root_files = my_alg.execute_on(scan_points + exmc_kk)