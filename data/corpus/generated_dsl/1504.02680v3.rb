# frozen_string_literal: true

### Dataset preparation ###

# ---------------------------------------------------------------------------
# Real data: continuum R-scan points (BOSS 713) plus the 3500-3671 MeV region
# (BOSS 704 psi(2S) scan points that fall inside that window).
# ---------------------------------------------------------------------------
data_2232 = DatasetManager.real_data.find("713_Rscan_2232")   # sqrt(s) ~ 2232.4 MeV
data_2400 = DatasetManager.real_data.find("713_Rscan_2396")   # sqrt(s) ~ 2396.4 MeV (closest R-scan point to 2400)
data_2800 = DatasetManager.real_data.find("713_Rscan_2800")   # sqrt(s) = 2800.0 MeV
data_3050 = DatasetManager.real_data.find("713_Rscan_3000")   # sqrt(s) = 3000.0 MeV (3050/3060 region)
data_3060 = DatasetManager.real_data.find("713_Rscan_3020")   # sqrt(s) = 3020.0 MeV (3050/3060 region)
data_3080 = DatasetManager.real_data.find("713_Rscan_3080")   # sqrt(s) = 3080.0 MeV
data_3580 = DatasetManager.real_data.find("704_psip_scan_1")  # sqrt(s) ~ 3581.5 MeV (3500-3671 MeV)
data_3670 = DatasetManager.real_data.find("704_psip_scan_2")  # sqrt(s) ~ 3670.2 MeV (3500-3671 MeV)

data_points = [data_2232, data_2400, data_2800, data_3050,
               data_3060, data_3080, data_3580, data_3670]

# Matching inclusive MC samples for every scan point
incMC_points = [
  DatasetManager.inclusive_mc.find("713_Rscan_2232"),
  DatasetManager.inclusive_mc.find("713_Rscan_2396"),
  DatasetManager.inclusive_mc.find("713_Rscan_2800"),
  DatasetManager.inclusive_mc.find("713_Rscan_3000"),
  DatasetManager.inclusive_mc.find("713_Rscan_3020"),
  DatasetManager.inclusive_mc.find("713_Rscan_3080"),
  DatasetManager.inclusive_mc.find("704_psip_scan_1"),
  DatasetManager.inclusive_mc.find("704_psip_scan_2")
]

# ---------------------------------------------------------------------------
# Signal decay card: ConExc generator (continuum / R-scan).
# ISR up to second order convoluted with the measured sigma0(m) is handled by
# ConExc; mode index 0 = p pbar.  The per-scan-point sqrt(s) ("Particle vpho")
# is injected automatically by the DSL for a multi-energy scan, so it is
# deliberately omitted here.
# ---------------------------------------------------------------------------
decay_card_signal_ppbar = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 0;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Background decay cards (KKMC + psi(4260) top-mother convention, since the
# description does not specify ISR modelling for the backgrounds).
# ---------------------------------------------------------------------------
decay_card_bkg_pipi = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_KpKm = <<~DECAYCARD
    Decay psi(4260)
    1.000 K+ K- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_ppbarpi0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 p+ anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_ppbarpi0pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 p+ anti-p- pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_LambdaLambdabar = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Signal MC: one 200k-event ConExc sample per scan energy point
# ---------------------------------------------------------------------------
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ppbar_signal_conexc"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_signal_ppbar
  config.cross_section = :default
end

# ---------------------------------------------------------------------------
# Background MC: 50k events per process per scan energy point
# ---------------------------------------------------------------------------
exMCs_bkg_pipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_pipi"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_pipi
  config.cross_section = :default
end

exMCs_bkg_KpKm = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_KpKm"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_KpKm
  config.cross_section = :default
end

exMCs_bkg_ppbarpi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_ppbarpi0"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_ppbarpi0
  config.cross_section = :default
end

exMCs_bkg_ppbarpi0pi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_ppbarpi0pi0"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_ppbarpi0pi0
  config.cross_section = :default
end

exMCs_bkg_LambdaLambdabar = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_LambdaLambdabar"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_LambdaLambdabar
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PpbarScan"
ppbar_alg = Algorithm.new(alg_name)
ppbar_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 2.232]})  # placeholder; per-scan-point beam energy injected at execution
         .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side procedures that cannot be expressed with the current DSL constructs
ppbar_alg
  .note(:ecms_per_scan_point, "ECMS varies across the R-scan (2232.4-3671.0 MeV); the constant declared here is a placeholder and the per-point beam energy is used at execution.")
  .note(:bhabha_veto, "E/p < 0.5 required for both proton and anti-proton candidates to suppress Bhabha (e+e- -> e+e-) contamination; applied after PID on the proton-candidate tracks (no dedicated DSL PID cut exists).")
  .note(:cosmic_veto, "TOF cosmic veto |T1 - T2| < 4 ns between the two proton-candidate tracks; applied before the kinematic fit.")
  .note(:proton_polar_angle_cut, "cos(theta_p) < 0.8 required for sqrt(s) > 2400 MeV; energy-dependent, applied per scan point.")
  .note(:ppbar_opening_angle_cut, "p pbar opening angle > 178 deg for sqrt(s) <= 2400 MeV and > 179 deg for sqrt(s) > 2400 MeV; energy-dependent, applied per scan point.")
  .note(:yield_extraction, "The published yield is obtained by event counting in the |p_meas - p_exp| < 5 sigma_p window without a kinematic fit; the 4C kinematic fit below is the BOSS-side endpoint only.")

# Charged-track, PID and kinematic-fit chain (common to all scan points)
event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==1"   # exactly one positively charged track
    nChrn     "==1"   # exactly one negatively charged track
    nNet      "==0"   # net charge zero
  }
  .pid(method: :probability) {   # PID combining dE/dx and TOF via the probability method
    prob_cut 0.001               # probability threshold 0.001
    identify :proton, against: [:kaon, :pion]  # one p+ and one anti-p- vs K and pi hypotheses
    nprp "==1"                   # one proton
    nprm "==1"                   # one anti-proton
  }
  # Nominal endpoint: 4C kinematic fit on the p pbar pair
  .kinematic_fit([:prp, :prm]) {
    nominal                    # nominal fit — corrected four-momenta are saved
    constrain_four_momentum    # 4C energy-momentum constraint to the CMS
    chi2_cut 200               # chi2 < 200 (loose BOSS cut; optimal cut made in ROOT)
  }

# Generate the algorithm for the p pbar signal process
ppbar_alg.with_decay_card(decay_card_signal_ppbar).apply(event_selection)

# Execute on all real data, inclusive MC and exclusive MC samples
root_files = ppbar_alg.execute_on(
  data_points + incMC_points +
  exMCs_signal +
  exMCs_bkg_pipi + exMCs_bkg_KpKm + exMCs_bkg_ppbarpi0 +
  exMCs_bkg_ppbarpi0pi0 + exMCs_bkg_LambdaLambdabar
)