### Dataset description ###
# ψ(3770) at √s = 3.773 GeV: real data and inclusive MC (BOSS 7.1.2 sample name 712_3773)
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process e+e- → D+D- (EvtGen format).
# Signal D+ → K+ π+ π- π0 ; tag-side D- in the three Cabibbo-favoured modes
# D- → K+π-π-, D- → K_S0π-, D- → K+π-π-π0 with relative fractions 0.340 / 0.270 / 0.390.
decay_card_dp_kpipipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 K+ pi+ pi- pi0  PHSP;
    Enddecay

    Decay D-
    0.3400 K+ pi- pi-       PHSP;
    0.2700 K_S0 pi-         PHSP;
    0.3900 K+ pi- pi- pi0   PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the D+ → K+ω, ω → π+π-π0 (Dalitz) signal sub-channel,
# using the same three D- tag modes.
decay_card_dp_komega = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 K+ omega  PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0  OMEGA_DALITZ;
    Enddecay

    Decay D-
    0.3400 K+ pi- pi-       PHSP;
    0.2700 K_S0 pi-         PHSP;
    0.3900 K+ pi- pi- pi0   PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the two signal modes
exMC_dp_kpipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKPiPiPi0"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_dp_kpipipi0
  config.cross_section   = :default
end

exMC_dp_komega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKOmega"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_dp_komega
  config.cross_section   = :default
end

### Tag-based event selection (double tag at ψ(3770)) ###
alg_name = "DpDmTagKPiPiPi0"
tag_alg  = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.773]})
       # tag mBC and ΔE are stored unconditionally (store-not-cut); the mode-dependent
       # ΔE windows are applied later in the ROOT mBC_tag vs mBC_sig 2D fit.
       .note(:tag_delta_e_window, "tag ΔE windows (-25, 25) MeV for tag modes without pi0 and
             (-55, 40) MeV for tag modes with pi0 are mode dependent and are applied in the ROOT
             mBC_tag vs mBC_sig 2D fit on the unconditionally stored tag ΔE/mBC, not at BOSS level")

# Tag side: D- from pre-stored DTag candidates, three Cabibbo-favoured modes, charge pinned to D-.
tag_alg.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0   # charge-conjugate of the canonical D+ mode names
  t.charm -1                                     # tag the D- side
end

# Signal side: D+ → K+π+π-π0 (incl. D+ → K+ω, ω → π+π-π0) built from the tag's remaining tracks/showers.
tag_alg.signal_side do |s|
  s.photons 2                       # exactly two good photons (π0 → γγ)
  s.charged(kp: 1, pip: 2, pim: 1)  # exactly one K+, two π+, one π- on the signal side
  s.require_charge 1                # net charge +1
  s.min_photon_angle 10.0           # photon-track isolation (DTagTool isGoodShower default)
end

# Kinematic fit: 4-momentum conservation plus γγ invariant mass constrained to the nominal π0 mass.
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200                    # loose χ² cut; tightened in the ROOT analysis
end

tag_alg.with_decay_card(decay_card_dp_kpipipi0).apply   # apply takes no Selection argument

root_files = tag_alg.execute_on([data_3773, incMC_3773, exMC_dp_kpipipi0, exMC_dp_komega])