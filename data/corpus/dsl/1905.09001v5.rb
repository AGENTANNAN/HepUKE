# BESIII DSL: Measurement of Proton Electromagnetic Form Factors in e+e- → ppbar
# ArXiv: 1905.09001v5
# 22 energy scan points from 2.00 to 3.08 GeV via ConExc mode 0 (p pbar)

### Dataset preparation ###

# ConExc decay card for e+e- → ppbar (mode 0 = p pbar, range 1.877-4.500 GeV)
# Omit Particle vpho for multi-energy scan — DSL injects per-point ECMS
decay_card_ppbar = <<~DECAYCARD
  Decay vpho
  1 ConExc 0;
  Enddecay
  Decay vhdr
  1 p+ anti-p- PHSP;
  Enddecay
  End
DECAYCARD

# 22 R-scan data points used in the analysis
ppbar_scan_points = [
  DatasetManager.real_data.find("713_Rscan_2000"),
  DatasetManager.real_data.find("713_Rscan_2050"),
  DatasetManager.real_data.find("713_Rscan_2100"),
  DatasetManager.real_data.find("713_Rscan_2125"),
  DatasetManager.real_data.find("713_Rscan_2150"),
  DatasetManager.real_data.find("713_Rscan_2175"),
  DatasetManager.real_data.find("713_Rscan_2200"),
  DatasetManager.real_data.find("713_Rscan_2232"),
  DatasetManager.real_data.find("713_Rscan_2309"),
  DatasetManager.real_data.find("713_Rscan_2386"),
  DatasetManager.real_data.find("713_Rscan_2396"),
  DatasetManager.real_data.find("713_Rscan_2500"),
  DatasetManager.real_data.find("713_Rscan_2644"),
  DatasetManager.real_data.find("713_Rscan_2646"),
  DatasetManager.real_data.find("713_Rscan_2700"),
  DatasetManager.real_data.find("713_Rscan_2800"),
  DatasetManager.real_data.find("713_Rscan_2900"),
  DatasetManager.real_data.find("713_Rscan_2950"),
  DatasetManager.real_data.find("713_Rscan_2981"),
  DatasetManager.real_data.find("713_Rscan_3000"),
  DatasetManager.real_data.find("713_Rscan_3020"),
  DatasetManager.real_data.find("713_Rscan_3080"),
]

# Signal exclusive MC via ConExc at each scan point
sig_ppbar = DatasetManager.create_exclusive_mc_for(ppbar_scan_points) do |config|
  config.sample_name   = "sig_conexc_ppbar"
  config.events        = 500_000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

### Event selection (BOSS) ###

alg_ppbar = Algorithm.new("PPbarFormFactor")
alg_ppbar.set_header(["PPbarFormFactorAlg/PPbarFormFactor.h"])
          .set_constant({ "ECMS" => [:double, 2.5] })

sel_ppbar = Selection.new
  # Exactly two charged tracks with opposite charge
  .select_track {
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
  }
  # Proton/antiproton identification via TOF and dE/dx (probability method)
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  # Remove identified proton/antiproton from charged track lists
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  # 4C kinematic fit constraining ppbar system to CMS energy
  .kinematic_fit([:prp, :prm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_ppbar
  .note(:vertex_fit_quality, "vertex fit of the two tracks to common vertex with chi2 < 100 before kinematic fit")
  .note(:opening_angle_cut, "opening angle between proton and antiproton in e+e- c.m. frame > 170 deg at 2.00/2.05 GeV, > 175 deg at 2.100-2.309 GeV, > 178 deg at 2.386-3.080 GeV")
  .note(:cosmic_rejection, "cosmic-ray background rejected by |T_trk1 - T_trk2| < 4 ns from TOF system")
  .note(:momentum_window, "tracks required within asymmetric momentum window (p_mean - 4sigma) < p < (p_mean + 3sigma); p_mean and sigma from fit to momentum distribution at each energy point")
  .note(:ep_cut, "E/p cut applied at energies > 2.150 GeV to suppress Bhabha events; E/p = E_EMC / p_MDC")
  .note(:iterative_mc_tuning, "iterative MC tuning (3 iterations): sigma_ppbar and |GE/GM| results fed back into CONEXC generation until stable within 1%")
  .note(:mc_model_uncertainty, "PHOKHARA-based alternative MC model used for systematic uncertainty estimation")
  .note(:angular_selection_range, "angular analysis restricted to |cos(theta)| < 0.8 due to TOF/EMC barrel-endcap gap efficiency")
  .note(:qed_background, "QED background from BABAYAGA and inclusive hadronic MC (CONEXC); estimated contamination < 0.5% and neglected")
  .with_decay_card(decay_card_ppbar)
  .apply(sel_ppbar)

alg_ppbar.execute_on(ppbar_scan_points + sig_ppbar)