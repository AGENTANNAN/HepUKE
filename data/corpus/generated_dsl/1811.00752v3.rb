### Dataset description ###
data_4178  = DatasetManager.real_data.find("703_4180")     # psi(4260) real data at sqrt(s) = 4.178 GeV (3.19 fb^-1)
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")  # corresponding inclusive MC sample

### Decay card for the signal process: psi(4260) -> D_s*+ D_s-, D_s*+ -> D_s+ gamma, D_s+ -> p nbar ###
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 D_s+ gamma PHSP;
    Enddecay

    Decay D_s+
    1.000 p+ anti-n- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 4 million signal events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4178_Dsstar_Ds_gamma_pnbar"
  config.related_dataset = data_4178
  config.events          = 4_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based analysis (TagAnalysis) ###
alg = TagAnalysis.new("DsStarDsPnbar")
alg.set_header(["DsStarDsPnbarAlg/DsStarDsPnbar.h"])
   .set_constant({"ECMS" => [:double, 4.178]})
   .set_alias({"std::vector<double>" => "Vdouble"})
   .with_decay_card(decay_card_signal)

# Tag side: reconstruct D_s- via the 11 single-tag modes
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,       # K_S K
          :DstoKKPi,      # K K pi
          :DstoKsKPi0,    # K_S K pi0
          :DstoKKPiPi0,   # K K pi pi0
          :DstoKsKPiPi,   # K_S K pi pi
          :DstoPiPiPi,    # pi pi pi
          :DstoPiEta,     # pi eta
          :DstoPiPi0Eta,  # pi pi0 eta
          :DstoPiPiPiPi,  # pi pi pi pi
          :DstoPiPiRhoGam,# pi pi rho gamma
          :DstoKPiPi      # K pi pi
  t.charm -1   # pin the tagged side to D_s- (charge -1)
end

# Signal side (what the tag did not use): D_s*+ -> D_s+ gamma, D_s+ -> p nbar (nbar missing)
alg.signal_side do |s|
  s.photons 1                     # one good photon (from D_s*+ -> D_s+ gamma)
  s.min_photon_angle 10.0         # photon angle to charged tracks > 10 degrees
  s.charged(prp: 1)               # exactly one charged track, charge +1 -> proton
  s.missing :n_bar                # the anti-neutron is the missing particle
end

# Kinematic fit: constrain D_s-, D_s+, D_s* masses and the initial four-momentum,
# with the anti-neutron treated as missing; best photon hypothesis selected by chi2 < 200
alg.fit do |f|
  f.constrain_four_momentum                                                      # initial-state 4-momentum constraint
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D_s-")               # tag side D_s- mass
  f.invariant_mass_of(:prp, :n_bar).constrain_to_nominal_mass_of(:"D_s+")        # D_s+ -> p nbar mass
  f.invariant_mass_of(:prp, :n_bar, :gamma).constrain_to_nominal_mass_of(:"D_s*+") # D_s*+ -> D_s+ gamma mass
  f.chi2_cut 200                                                                 # chi2 < 200
end

alg.apply   # validate + render the tag spec (no Selection argument)
alg.execute_on([data_4178, incMC_4178, exMC_signal])