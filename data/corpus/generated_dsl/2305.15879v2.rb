### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data at 3.773 GeV (~2.93 fb^-1)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC sample

# Decay card for the signal process: psi(3770) -> D+ D-, D+ -> K_S0 pi+ pi0 pi0
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 pi+ pi0 pi0 PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for the D+ -> K_S0 pi+ pi0 pi0 signal mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dp_kspi0pi0"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (double tag: D- tag + D+ signal) ###
alg_name = "DpKsPiPi0Pi0"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: six hadronic D- tag modes
alg.tag_side(:Dm) do |t|
  t.modes :DmtoKPiPi, :DmtoKPiPiPi0, :DmtoKsPi, :DmtoKsPiPi0, :DmtoKKPi, :DmtoKPiPiPi
  t.charm -1                                    # tag the D- side
  t.window :deltaE, min: -0.025, max: 0.020     # tag-side DeltaE window [-25, 20] MeV
end

# Signal side: D+ -> K_S0 pi+ pi0 pi0 (K_S0 -> pi+ pi-, two pi0 -> four photons)
alg.signal_side do |s|
  s.charged(pip: 2, pim: 1)   # K_S0 -> pi+ pi- plus the direct pi+; net charge +1
  s.photons 4                 # four photons from the two pi0
  s.require_charge(1)
end

# Kinematic fit: 4C four-momentum conservation + M(gamma gamma) = m(pi0), chi2 < 200
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures preserved for downstream handling
alg.note(:background_veto, "K_S0 veto: reject events with M(pi0 pi0) in [0.46, 0.52] GeV/c^2")
   .note(:background_veto, "D0 D0bar wrong-combination veto: tag M_BC required outside [1.862, 1.870] GeV")
   .note(:deltaE_window, "signal-side DeltaE window [-40, 20] MeV applied in ROOT")
   .note(:fit_amplitude, "for the amplitude analysis restricted to the D- -> K+ pi- pi- tag, a 4C fit constraining the D+, K_S0 and pi0 masses is used")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])