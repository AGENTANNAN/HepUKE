### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # 20.3 fb^-1 psi(3770) data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # corresponding inclusive MC

# --- decay card: psi(3770) -> D+ D-, D+ -> pi+ pi+ pi- eta, eta -> gamma gamma, D- -> K+ pi- pi- ---
decay_card_modeI = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 pi+ pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# --- decay card: psi(3770) -> D+ D-, D+ -> pi+ pi0 pi0 eta, pi0/eta -> gamma gamma, D- -> K+ pi- pi- ---
decay_card_modeII = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 pi+ pi0 pi0 eta PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# --- exclusive MC: 500k events for each signal mode ---
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dpto3piEta"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dptopipi0pi0Eta"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (TagAnalysis: D- tag on one side, D+ signal on the other) ###

# ---------- Mode I: D+ -> pi+ pi+ pi- eta, eta -> gamma gamma ----------
alg_modeI = TagAnalysis.new("DpTo3PiEtaDT")
alg_modeI.set_header(["DpTo3PiEtaDTAlg/DpTo3PiEtaDT.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .note(:background_veto, "K_S veto: require M(pi+pi-) outside [0.478, 0.517] GeV/c2; eta' veto: require M(pi+pi-eta) outside [0.8, 1.0] GeV/c2.")
         .note(:candidate_selection, "when several candidates survive, choose the one with minimum |dE_sig| and, for this channel, minimum |dE_tag|.")

# Tag the D- in the five hadronic modes (charm -1 selects the D- side)
alg_modeI.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: the tracks/showers left over opposite the tag
alg_modeI.signal_side do |s|
  s.charged(pip: 2, pim: 1)   # 2 pi+ and 1 pi- on the signal side
  s.require_charge 1          # net charge +1
  s.photons 2                 # two photons from eta -> gamma gamma
  s.min_photon_angle 10.0     # photon angle to charged tracks > 10 deg
  s.min_photon_energy 0.025   # photon energy > 25 MeV
end

# Fit: 4-momentum conservation + gamma-gamma mass constrained to the eta
alg_modeI.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 100
end

alg_modeI.with_decay_card(decay_card_modeI).apply
alg_modeI.execute_on([data_3773, incMC_3773, exMC_modeI])

# ---------- Mode II: D+ -> pi+ pi0 pi0 eta, pi0/eta -> gamma gamma ----------
alg_modeII = TagAnalysis.new("DpToPiPi0Pi0EtaDT")
alg_modeII.set_header(["DpToPiPi0Pi0EtaDTAlg/DpToPiPi0Pi0EtaDT.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:background_veto, "K_S veto: require M(pi0pi0) outside [0.448, 0.534] GeV/c2; eta' veto on pi0 eta; D+ -> pi+pi0pi0pi0 background rejected with dE in (-0.1, 0.1) GeV and M_BC in [1.83, 1.89] GeV/c2.")
          .note(:candidate_selection, "when several candidates survive, choose the one with minimum |dE_sig|.")

# Same shared tag side
alg_modeII.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: one pi+ and six photons (pi0 -> gg, pi0 -> gg, eta -> gg)
alg_modeII.signal_side do |s|
  s.charged(pip: 1)           # 1 pi+ on the signal side
  s.require_charge 1          # net charge +1
  s.photons 6                 # six photons from pi0 pi0 eta -> gamma gamma
  s.min_photon_angle 10.0     # photon angle to charged tracks > 10 deg
  s.min_photon_energy 0.025   # photon energy > 25 MeV
end

# Fit: 4-momentum conservation + two gamma-gamma constraints to pi0 + one to eta
alg_modeII.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 50
end

alg_modeII.with_decay_card(decay_card_modeII).apply
alg_modeII.execute_on([data_3773, incMC_3773, exMC_modeII])