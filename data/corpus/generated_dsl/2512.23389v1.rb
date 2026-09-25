### Dataset preparation ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data (~20.3 fb⁻¹)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# Exclusive signal MC: ψ(3770) → D0 D0bar, D0 → K_S0 π0 η, D0bar → K+ π- π0
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 pi0 eta PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0toKsPi0Eta_D0bartoKPiPi0"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (tag-based, D0 D0bar double tag) ###
alg = TagAnalysis.new("D0D0barDTag")
alg.set_header(["D0D0barDTagAlg/D0D0barDTag.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)

# Tag side: D0bar (charm -1) reconstructed in the three hadronic modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :mBC, min: 1.859, max: 1.873   # tag M_BC selection (explicitly requested)
end

# Signal side: D0 → K_S0 π0 η with K_S0 → π+π-, π0 → γγ, η → γγ
alg.signal_side do |s|
  s.photons 4..4                 # the four photons from π0 and η, no additional photons
  s.charged(pip: 1, pim: 1)      # π+π- from K_S0
end

# 4C kinematic fit: four-momentum conservation plus nominal D0bar, D0, K_S0, π0, η masses
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)                # tag D0bar mass
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)         # K_S0 → π+π-
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)      # π0 → γγ
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)      # η → γγ
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma, :gamma, :gamma)
   .constrain_to_nominal_mass_of(:D0)                                         # signal D0 mass
  f.chi2_cut 200
  f.store_fitted_momenta                                                      # store fitted four-momenta
end

# BOSS-side procedures that have no DSL construct
alg.note(:background_veto,
         "Veto events in which any pair of the four π0/η photons has an invariant mass in 0.115-0.150 GeV/c², to suppress π0-η photon mis-pairing")
alg.note(:signal_side_selection,
         "Signal-side D0 selection: 1.858 < M_BC < 1.871 GeV/c² and -50 < ΔE < 45 MeV; among multiple signal candidates the one with minimum |ΔE| is chosen (stored mBC/ΔE are windowed and the best candidate is picked in the ROOT analysis)")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])