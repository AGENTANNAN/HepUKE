# Study of eta -> l+ l- (l = e, mu) via J/psi -> gamma eta'
#   [arXiv:2512.07144]
#
# Decay chain: J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> l+ l-
# Using (10087 +/- 44) x 10^6 J/psi events at sqrt(s) = 3.097 GeV.
# Two signal modes: (I) eta -> mu+ mu-, (II) eta -> e+ e-.
# This is an ordinary (non-tag) analysis: Selection built from tracks/photons/PID.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Inclusive J/psi MC

# ---------------------------------------------------------------- decay cards
# Sub-decay for eta -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'  PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta  PHSP;
  Enddecay

  Decay eta
  1.0000 mu+ mu-  PHSP;
  Enddecay

  End
DECAYCARD

# Sub-decay for eta -> e+ e-
decay_card_ee = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'  PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta  PHSP;
  Enddecay

  Decay eta
  1.0000 e+ e-  PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for mode I: eta -> mu+ mu-
exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_gamma_etap_pipi_eta_mumu"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

# Exclusive MC for mode II: eta -> e+ e-
exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_gamma_etap_pipi_eta_ee"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

### Event selection (BOSS) — Mode I: eta -> mu+ mu- ###
alg_mumu = Algorithm.new("JpsiEtaP2PiPiEta2MuMu")
alg_mumu.set_header(["JpsiEtaP2PiPiEta2MuMuAlg/JpsiEtaP2PiPiEta2MuMu.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_mumu = Selection.new
sel_mumu.select_track {
           cos_theta 0.93          # |cos(theta)| < 0.93
           Vz        100.0        # |Vz| < 10 cm
           Vr        10.0          # |Vr| < 1 cm
           nChrp     ">=2"
           nChrn     ">=2"
           nNet      "==0"
         }
        .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0
           energyThreshold_b 0.025    # >=25 MeV in barrel
           energyThreshold_e 0.050    # >=50 MeV in endcap
           nGam              ">=1"
         }
        .pid(method: :probability) {
           prob_cut 0.001
           identify :muon, against: [:pion, :kaon]
           identify :pion, against: [:kaon]
           nmup     ">=1"
           nmum     ">=1"
         }
        .remove([:mup <= :chrgp])
        .remove([:mum <= :chrgn])
        .assign({:chrgp => :pip, :chrgn => :pim})
        .kinematic_fit([:gamma, :pip, :pip, :pim, :pim]) {
           nominal
           constrain_four_momentum
           chi2_cut 40            # chi2_4C+PID < 40
         }

alg_mumu
  .note(:pid_chi2_sum,
        "chi2_4C+PID = chi2_4C + chi2_PID used for best-candidate selection. " \
        "The PID chi2 sum is over the four charged track candidates; " \
        "the hypothesis with smallest chi2_4C+PID is retained.")
  .note(:background_veto,
        "Competing hypothesis: 4C fit under J/psi -> gamma 2(pi+ pi-) " \
        "(i.e., J/psi -> gamma pi+ pi- pi+ pi-). " \
        "Require chi2_4C+PID(pi+ pi- pi+ pi-) > chi2_4C+PID(pi+ pi- mu+ mu-) " \
        "to reject J/psi -> gamma 2(pi+ pi-) background. " \
        "Muon counter information is NOT used due to low muon momentum.")
  .note(:tagging_chain,
        "Analysis uses the novel approach of tagging eta' -> pi+ pi- eta in J/psi radiative decays. " \
        "M(pi+ pi- mu+ mu-) in [0.945, 0.970] GeV/c^2 selects eta' mass window. " \
        "The eta -> mu+ mu- signal yield is extracted from an extended unbinned ML fit " \
        "to M(mu+ mu-), with signal shape from MC and backgrounds from " \
        "eta' -> pi+ pi- pi+ pi- (fixed 0.8 +/- 0.1), eta' -> pi+ pi- mu+ mu- (fixed 15.0 +/- 3.0), " \
        "and J/psi -> gamma pi+ pi- pi+ pi- (free).")
  .note(:eta_ee_upper_limit,
        "For eta -> e+ e-: No events observed in the signal region. " \
        "Signal region: M(e+e-) in [0.534, 0.560] GeV/c^2 (3.5 sigma mass window). " \
        "Sideband regions: [0.508, 0.534] and [0.560, 0.586] GeV/c^2. " \
        "N_sig=0, N_bkg=1.5 estimated from sidebands. Upper limit at 90% CL: " \
        "B(eta -> e+ e-) < 2.2 x 10^-7. Detection efficiency (27.48 +/- 0.04)%.")
  .note(:efficiency_calculation,
        "BF(eta -> l+ l-) = N_obs / (N_J/psi * B_int * epsilon). " \
        "B_int = B(J/psi -> gamma eta') * B(eta' -> pi+ pi- eta) from PDG. " \
        "Detection efficiency from signal MC: (29.19 +/- 0.07)% for mu+ mu-, " \
        "(27.48 +/- 0.04)% for e+ e- (statistical only).")
  .note(:pion_tracking_sys,
        "Pion tracking efficiency correction uses J/psi -> pi+ pi- pi0 control sample. " \
        "Muon tracks weighted by pion tracking corrections since low-momentum muons " \
        "behave similarly to pions in the MDC.")
  .note(:electron_tracking_sys,
        "Electron tracking efficiency correction uses e+ e- -> e+ e- gamma and " \
        "J/psi -> e+ e- gamma_FSR control samples.")
  .with_decay_card(decay_card_mumu)
  .apply(sel_mumu)

### Event selection (BOSS) — Mode II: eta -> e+ e- ###
alg_ee = Algorithm.new("JpsiEtaP2PiPiEta2EE")
alg_ee.set_header(["JpsiEtaP2PiPiEta2EEAlg/JpsiEtaP2PiPiEta2EE.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

sel_ee = Selection.new
sel_ee.select_track {
         cos_theta 0.93
         Vz        100.0
         Vr        10.0
         nChrp     ">=2"
         nChrn     ">=2"
         nNet      "==0"
       }
      .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         angle_to_track    10.0
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         nGam              ">=1"
       }
      .pid(method: :probability) {
         prob_cut 0.001
         identify :electron, against: [:pion, :kaon]
         identify :pion, against: [:kaon]
         nep      ">=1"
         nem      ">=1"
       }
      .remove([:ep <= :chrgp])
      .remove([:em <= :chrgn])
      .assign({:chrgp => :pip, :chrgn => :pim})
      .kinematic_fit([:gamma, :pip, :pip, :pim, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 40
       }

alg_ee
  .note(:pid_chi2_sum,
        "chi2_4C+PID = chi2_4C + chi2_PID used for best-candidate selection. " \
        "The PID chi2 sum is over the four charged track candidates; " \
        "the hypothesis with smallest chi2_4C+PID is retained.")
  .note(:eta_ee_upper_limit,
        "No events observed in the signal region. " \
        "Signal region: M(e+e-) in [0.534, 0.560] GeV/c^2 (3.5 sigma mass window). " \
        "Background from J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma e+ e- " \
        "(0.5 +/- 0.3) and eta' -> pi+ pi- e+ e- (1.0 +/- 0.2) from dedicated MC. " \
        "Upper limit at 90% CL: B(eta -> e+ e-) < 2.2 x 10^-7.")
  .with_decay_card(decay_card_ee)
  .apply(sel_ee)

alg_mumu.execute_on([jpsi_data, jpsi_incMC, exMC_mumu])
alg_ee.execute_on([jpsi_data, jpsi_incMC, exMC_ee])