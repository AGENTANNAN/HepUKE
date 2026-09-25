### Dataset description ###
# J/psi(3097) real data and the corresponding inclusive MC (sample name 708_3097)
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# --- Signal decay cards (EvtGen format) ---
# Mode I: J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Mode II: J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC: 200k events for each eta decay mode ---
exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_etap_gamma_pipim_eta_mumu"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_etap_gamma_pipim_eta_ee"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ===================== Mode I: eta -> mu+ mu- =====================
alg_mumu = Algorithm.new("JpsiGammaEtaPrimeEtaToMuMu")
alg_mumu.set_header(["JpsiGammaEtaPrimeEtaToMuMuAlg/JpsiGammaEtaPrimeEtaToMuMu.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_mumu = Selection.new
  .select_track {                 # charged track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     ">=2"               # at least 2 positive tracks
    nChrn     ">=2"               # at least 2 negative tracks
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # photon selection
    tdc_emc_start     0           # EMC TDC start
    tdc_emc_end       14          # EMC TDC end
    angle_to_track    10.0        # angle to nearest charged track > 10 degrees
    energyThreshold_b 0.025       # barrel energy threshold 25 MeV
    energyThreshold_e 0.050       # endcap energy threshold 50 MeV
    nGam              ">=1"       # at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    # muon identification (e/mu are not separable by identify(); use the lepton finder)
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]   # pi / K separation
    nlp  ">=1"                    # at least one mu+
    nlm  ">=1"                    # at least one mu-
    npip ">=1"                    # at least one pi+
    npim ">=1"                    # at least one pi-
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])    # remove identified leptons from the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})  # remaining tracks -> pi+ / pi-
  # competing gamma pi+ pi- mu+ mu- hypothesis: stores chi2 (no nominal, no chi2 cut)
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
    constrain_four_momentum
  }
  # nominal 4C fit under the all-pion hypothesis: every charged track treated as a pion
  .assign({:lp => :pip, :lm => :pim})        # fold the lepton tracks back in as pions
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }

alg_mumu.note(:best_candidate_selection,
              "best candidate combination chosen by the sum of the 4C kinematic-fit chi2 and the PID chi2")
        .note(:background_veto,
              "in the muon mode the all-pion gamma pi+ pi- pi+ pi- chi2 must exceed the gamma pi+ pi- mu+ mu- chi2 in order to reject J/psi -> gamma 2(pi+ pi-) background; the competing chi2 is stored by the non-nominal fit")
        .with_decay_card(decay_card_mumu).apply(sel_mumu)

alg_mumu.execute_on([jpsi_data, jpsi_incMC, exMC_mumu])

# ===================== Mode II: eta -> e+ e- =====================
alg_ee = Algorithm.new("JpsiGammaEtaPrimeEtaToEE")
alg_ee.set_header(["JpsiGammaEtaPrimeEtaToEEAlg/JpsiGammaEtaPrimeEtaToEE.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_ee = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
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
    # electron identification (e/mu are not separable by identify(); use the lepton finder)
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]   # pi / K separation
    nlp  ">=1"                    # at least one e+
    nlm  ">=1"                    # at least one e-
    npip ">=1"
    npim ">=1"
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  # nominal 4C fit under the all-pion hypothesis: every charged track treated as a pion
  .assign({:lp => :pip, :lm => :pim})
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }

alg_ee.note(:best_candidate_selection,
            "best candidate combination chosen by the sum of the 4C kinematic-fit chi2 and the PID chi2")
      .with_decay_card(decay_card_ee).apply(sel_ee)

alg_ee.execute_on([jpsi_data, jpsi_incMC, exMC_ee])