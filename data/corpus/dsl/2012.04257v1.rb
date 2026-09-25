# BESIII DSL: Observation of eta' -> pi+ pi- mu+ mu-
# Paper: 2012.04257v1
# J/psi -> gamma eta' at sqrt(s) = 3.097 GeV

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Decay cards
# ============================================================

# Signal MC: J/psi -> gamma eta', eta' -> pi+ pi- mu+ mu-
decay_card_sig = <<~DECAYCARD
  Decay J/psi
    1 gamma eta_prime VSS;
  Enddecay
  Decay eta_prime
    1 pi+ pi- mu+ mu- PHSP;
  Enddecay
  End
DECAYCARD

# Signal MC: J/psi -> gamma eta', eta' -> pi+ pi- pi+ pi- (competing background)
decay_card_4pi = <<~DECAYCARD
  Decay J/psi
    1 gamma eta_prime VSS;
  Enddecay
  Decay eta_prime
    1 pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC samples
sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_etap_pipimumu"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_sig
  config.cross_section   = :default
end

sig_4pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_etap_4pi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_4pi
  config.cross_section   = :default
end

# ============================================================
# Algorithm: eta' -> pi+ pi- mu+ mu- via J/psi -> gamma eta'
# 4 charged tracks (pi+, pi-, mu+, mu-) + at least 1 photon
# 4C kinematic fit under gamma pi+ pi- mu+ mu- with chi2 < 30
# Competing pi+pi-pi+pi- veto via chi2 comparison
# eta -> mu+mu- veto via |M(mu+mu-) - M_eta| > 0.02 GeV
# ============================================================

alg_etap = Algorithm.new("EtapPipimumu", version: '00-00-01')
alg_etap.set_header(["EtapPipimumuAlg/EtapPipimumu.h"])
         .set_constant({ "ECMS" => [:double, 3.097] })

etap_selection = Selection.new
  .select_track do
    nChrp ">=2"      # pi+, mu+
    nChrn ">=2"      # pi-, mu-
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    15.0     # > 15 degrees from nearest charged track
    nGam              ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]    # pion PID for pi+ and pi-
  end
  .pid do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.3,
                                   treat_as_muon_if_momentum_above: 0.3,
                                   treat_as_muon_if_energy_range: [0.05, 0.35]
    nlp "==0"; nlm ">=2"          # no electrons, at least two muons
  end
  # Nominal 4C kinematic fit: gamma pi+ pi- mu+ mu-
  .kinematic_fit([:gamma, :pip, :pim, :mup, :mum]) do
    constrain_four_momentum
    chi2_cut 30
    nominal
  end
  # Competing hypothesis: gamma pi+ pi- pi+ pi- (non-nominal, no chi2 cut)
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    # chi2 stored but no cut; competing chi2 comparison done in ROOT
  end

# Notes for inexpressible BOSS-side procedures
alg_etap.note(:combined_chi2, "Combined chi2 = chi2_4C + chi2_PID used for best candidate selection in ROOT")
alg_etap.note(:etap_mass_window, "eta' mass window (0.94, 0.97) GeV/c^2 applied in ROOT")
alg_etap.note(:eta_mumu_veto, "|M(mu+mu-) - M_eta| > 0.02 GeV/c^2 applied in ROOT (eta -> mu+mu- veto)")
alg_etap.note(:four_pi_veto, "chi2(gamma 4pi) < chi2(gamma pi+pi- mu+mu-) event discarded in ROOT")
alg_etap.note(:photon_from_etap, "Radiative photon from J/psi decay distinguished from eta' decay products by kinematic fit in ROOT")

alg_etap.with_decay_card(decay_card_sig).apply(etap_selection)
alg_etap.execute_on([jpsi_data, jpsi_incMC, sig_mc, sig_4pi])