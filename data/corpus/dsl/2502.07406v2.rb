# BESIII DSL: Search for e+e- -> K_S0 K_S0 h_c
# Paper: 2502.07406v2
# 13 energy points: 4.600 - 4.951 GeV
# Partial reconstruction: h_c -> gamma eta_c, eta_c detected only as recoil

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 13 energy point datasets (BOSS version / sample name)
# Paper labels vs. dataset entries: 4600->703_4600, 4612->706_4610, 4628->706_4620,
# 4641->706_4640, 4661->706_4660, 4682->706_4680, 4699->706_4700,
# 4740->707_4740, 4750->707_4750, 4781->707_4780, 4843->707_4840,
# 4918->707_4914, 4951->707_4946

ALL_ENERGY_KEYS = %w[703_4600 706_4610 706_4620 706_4640 706_4660
                      706_4680 706_4700 707_4740 707_4750 707_4780
                      707_4840 707_4914 707_4946]

all_data = ALL_ENERGY_KEYS.map { |k| DatasetManager.real_data.find(k) }
all_inc_mc = ALL_ENERGY_KEYS.map { |k| DatasetManager.inclusive_mc.find(k) }

# ============================================================
# Decay card: psi(4260) -> K_S0 K_S0 h_c, h_c -> gamma eta_c
# ============================================================

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 K_S0 K_S0 h_c PHSP;
  Enddecay
  Decay h_c
  1.000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.000 gamma gamma PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive signal MC (13 data points -> 13 MC samples)
sig_mcs = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_ksks_hc"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# ============================================================
# Algorithm: e+e- -> K_S0 K_S0 h_c (partial reconstruction)
#
# Reconstruct: K_S0(1) + K_S0(2) + gamma(4) from h_c decay
# Miss: eta_c(5) — detected as recoil mass peak
# ============================================================

algorithm = Algorithm.new("KsKshcPartialRec")

algorithm.set_header(["KsKshcPartialRec/KsKshcPartialRec.h"])

# ECMS varies across 13 energy points; representative value set here
# per-data-point fits use per-run MeasuredEcmsSvc
algorithm.set_constant({ "ECMS" => [:double, 4.600] })

algorithm.with_decay_card(decay_card)

event_selection = Selection.new
  # Charged tracks: two pi+, two pi- (for two K_S0 -> pi+ pi- decays)
  .select_track do
    nChrp "==2"
    nChrn "==2"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  # Reconstruct K_S0 candidates from oppositely charged pion pairs
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
  end
  # Photon selection: at least one good photon
  .select_photon do
    nGam ">=1"
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  end
  # Partial reconstruction: miss eta_c (recID 5), reconstruct K_S0(1), K_S0(2), gamma(4)
  # DecayCardResolver rec_id_list:
  #   0 => psi(4260)  (skip — replaced by P4_cms)
  #   1 => K_S0       (first)
  #   2 => K_S0       (second)
  #   3 => h_c
  #   4 => gamma      (E1 photon from h_c -> gamma eta_c)
  #   5 => eta_c      (missed — recoil)
  .partial_rec([1, 2, 3, 4]) do
    require_recoil_mass 2.94, 3.06
  end

algorithm.apply(event_selection)

algorithm.note(:ks_mass_window, "M(pi+ pi-) within m_K_S0 +/- 6 MeV/c^2")
algorithm.note(:ks_decay_length, "K_S0 decay length > 2 * vertex resolution")
algorithm.note(:photon_timing, "EMC time within [0, 700] ns of event start")
algorithm.note(:multi_ks_pair, "All K_S0 K_S0 combinations retained; each pion used once")
algorithm.note(:e1_photon_selection,
  "Photon with M_rec(gamma K_S0 K_S0) closest to m_eta_c selected as E1 photon")
algorithm.note(:multi_energy_fits,
  "13 energy points individually fitted with Umiss/M_rec; signal shape from MC " \
  "convolved with Gaussian (resolution from K_S0 K_S0 J/psi control sample); " \
  "background = 2nd-order Chebyshev")
algorithm.note(:control_sample,
  "e+e- -> K_S0 K_S0 J/psi control sample at 4.190-4.290 GeV used for " \
  "resolution calibration")

algorithm.execute_on(all_data + all_inc_mc + sig_mcs)