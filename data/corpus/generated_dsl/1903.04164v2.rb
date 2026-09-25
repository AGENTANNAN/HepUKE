# Core DSL classes and dependencies are loaded automatically at execution.
# Tag-based (D_s tag) analyses: e+e- -> D_s*± D_s∓ at sqrt(s) = 4.178 GeV (3.19 fb^-1).
# Two independent signal channels share the same tag side (single-tag D_s- hadronic modes):
#   (1) D_s+ -> K_S0 K+   (double tag; K_S0 -> pi+ pi- from the tag-unused tracks)
#   (2) D_s+ -> K_L0 K+   (single tag + missing K_L0; 4C kinematic fit)

### Dataset description ###
data_4178  = DatasetManager.real_data.find("703_4180")     # 4.178 GeV real data (3.19 fb^-1)
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")  # matching inclusive MC

### Decay cards (ConExc generator: e+e- -> D_s*± D_s∓, ISR modelled by ConExc) ###
# NOTE: the ConExc channel index for D_s*± D_s∓ is taken from the ConExc mode table;
# `Particle vpho <ECMS>` is omitted (injected automatically per energy point).
decay_card_ksk = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 20;
    Enddecay

    Decay D_s*+
    1.0 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.0 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.0 K_S0 K+ PHSP;
    Enddecay

    Decay D_s-
    1.0 K+ K- pi- PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_klk = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 20;
    Enddecay

    Decay D_s*+
    1.0 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.0 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.0 K_L0 K+ PHSP;
    Enddecay

    Decay D_s-
    1.0 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (ConExc, one per signal mode) ###
exMC_ksk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4178_ds_ks0k"
  config.related_dataset = data_4178
  config.events          = 100000
  config.decay_card      = decay_card_ksk
  config.cross_section   = :default
end

exMC_klk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4178_ds_kl0k"
  config.related_dataset = data_4178
  config.events          = 100000
  config.decay_card      = decay_card_klk
  config.cross_section   = :default
end

### ------------------------------------------------------------------ ###
### Channel 1: D_s+ -> K_S0 K+  (D_s- single tag + D_s+ -> K_S0 K+)   ###
### ------------------------------------------------------------------ ###
alg_ksk = TagAnalysis.new("DsToKsK")
alg_ksk.set_header(["DsToKsKAlg/DsToKsK.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
alg_ksk.with_decay_card(decay_card_ksk)

# Tag side: single-tag D_s- hadronic modes (the D_s- -> K_S0 K- mode is NOT declared,
# to avoid double counting against the D_s+ -> K_S0 K+ measurement).
alg_ksk.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoKsKPi0, :DstoKsKsPi,
          :DstoKKPiPiPi, :DstoKsKPiPi, :DstopiEta, :DstopiEtaPrime,
          :DstopiPiPiPi, :DstopiPi0, :DstoKKPiPi0, :DstoKsKPiPi0,
          :DstoKsKsPiPi
  t.charm -1                       # tag the D_s- side
end

# Signal side: one K+ (charge +1) plus the K_S0 -> pi+ pi- built from the tag-unused tracks.
alg_ksk.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)
  s.require_charge 1
end

# No constrained kinematic fit for this channel: only the K_S0 mass window is declared.
alg_ksk.fit do |f|
  f.invariant_mass_of(:pip, :pim).between(0.487, 0.511)   # K_S0 mass window (GeV/c^2)
end

alg_ksk.note(:tag_mode_exclusion,
             "the D_s- -> K_S0 K- single-tag mode is deliberately excluded from the tag "
             "mode list to avoid double counting against the D_s+ -> K_S0 K+ signal channel")
     .note(:k_s0_secondary_vertex,
             "the K_S0 -> pi+ pi- is reconstructed from the tag-unused tracks via a "
             "secondary-vertex fit and required in the mass window 0.487-0.511 GeV/c^2")
     .note(:extra_track_veto,
             "events with extra charged tracks satisfying |cos(theta)| < 0.93 and "
             "|Vz| < 20 cm are vetoed")
     .note(:no_kinematic_fit,
             "no constrained kinematic fit is applied for the D_s+ -> K_S0 K+ channel; the "
             "yield is extracted in ROOT from a 2D unbinned maximum-likelihood fit of "
             "M(K_S0 K+) versus M_tag")

alg_ksk.apply                       # tag analyses: apply takes NO Selection argument

### ------------------------------------------------------------------ ###
### Channel 2: D_s+ -> K_L0 K+  (D_s- single tag + missing K_L0)      ###
### ------------------------------------------------------------------ ###
alg_klk = TagAnalysis.new("DsToKlK")
alg_klk.set_header(["DsToKlKAlg/DsToKlK.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
alg_klk.with_decay_card(decay_card_klk)

# Identical tag side as channel 1 (single-tag D_s- hadronic modes).
alg_klk.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoKsKPi0, :DstoKsKsPi,
          :DstoKKPiPiPi, :DstoKsKPiPi, :DstopiEta, :DstopiEtaPrime,
          :DstopiPiPiPi, :DstopiPi0, :DstoKKPiPi0, :DstoKsKPiPi0,
          :DstoKsKsPiPi
  t.charm -1
end

# Signal side: one K+ (charge +1) plus the gamma from D_s* -> gamma D_s; K_L0 is missing.
alg_klk.signal_side do |s|
  s.photons 1
  s.charged(kp: 1)
  s.require_charge 1
  s.missing :K_L0                    # massive missing particle (K_L0)
end

# Nominal 4C fit: 4-momentum conservation + masses of the ST D_s-, the signal D_s+ and the
# D_s* intermediate state (signal-side D_s* -> gamma D_s hypothesis).
alg_klk.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D_s-")               # ST D_s-
  f.invariant_mass_of(:kp, :K_L0).constrain_to_nominal_mass_of(:"D_s+")          # signal D_s+
  f.invariant_mass_of(:kp, :K_L0, :gamma).constrain_to_nominal_mass_of(:"D_s*+") # signal-side D_s*
  f.chi2_cut 40
end

# Competing photon hypothesis: the gamma originates from the tag-side D_s*- -> gamma D_s-.
alg_klk.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D_s-")
  f.invariant_mass_of(:kp, :K_L0).constrain_to_nominal_mass_of(:"D_s+")
  f.invariant_mass_of(:tag1, :gamma).constrain_to_nominal_mass_of(:"D_s*-")      # tag-side D_s*
  f.chi2_cut 40
end

alg_klk.note(:ds_star_photon_choice,
             "the photon from D_s* -> gamma D_s is chosen by testing both the tag-side and "
             "signal-side D_s* hypotheses and keeping the smaller chi^2 (the two fit "
             "hypotheses above store both chi^2 values for the ROOT-level comparison)")
     .note(:extra_photon_veto,
             "extra photons with E > 250 MeV and an opening angle > 15 degrees to the missing "
             "particle are rejected")
     .note(:missing_mass_extraction,
             "the K_L0 signal is extracted from the missing-mass-squared "
             "MM^2 = (P_{e+e-} - P_{D_s-} - P_gamma - P_{K+})^2")
     .note(:peaking_background,
             "peaking backgrounds from D_s+ -> K_S0 K+ and D_s+ -> eta K+ are accounted for "
             "in the extraction")

alg_klk.apply

### Execute both algorithms on real data, inclusive MC and both signal MC samples ###
root_files_ksk = alg_ksk.execute_on([data_4178, incMC_4178, exMC_ksk, exMC_klk])
root_files_klk = alg_klk.execute_on([data_4178, incMC_4178, exMC_ksk, exMC_klk])