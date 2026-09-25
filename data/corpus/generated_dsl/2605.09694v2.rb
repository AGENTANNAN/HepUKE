### Dataset description ###
ds_data  = DatasetManager.real_data.find("703_4180")        # 4.178 GeV D_s* D_s data point (~7.33 fb^-1, 4.128-4.226 GeV)
ds_incMC = DatasetManager.inclusive_mc.find("703_4180")     # corresponding inclusive MC sample

# Decay card — signal mode I: D_s+ -> K_S0 K_S0 pi+ pi0 (K_S0 -> pi+pi-, pi0 -> gamma gamma)
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.0000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s+
    1.0000 K_S0 K_S0 pi+ pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — signal mode II: D_s+ -> K_S0 K+ pi0 pi0 (K_S0 -> pi+pi-, pi0 -> gamma gamma)
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.0000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s+
    1.0000 K_S0 K+ pi0 pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each of the two signal modes
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Ds_doubletag_KSKSpi_pi0_ModeI"
  config.related_dataset = ds_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Ds_doubletag_KSKpi_pi0pi0_ModeII"
  config.related_dataset = ds_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based double-tag analysis ###
# ---------- Mode I: D_s+ -> K_S0 K_S0 pi+ pi0 ----------
alg_modeI = TagAnalysis.new("DsTagModeI")
alg_modeI.set_header(["DsTagModeIAlg/DsTagModeI.h"])
         .set_constant({"ECMS" => [:double, 4.178]})
         .with_decay_card(decay_card_modeI)

# Tag side: the 16 hadronic D_s- modes (charm = -1)
alg_modeI.tag_side(:Ds) do |t|
  t.mode_group(:hadronic)   # all 16 hadronic D_s- tag modes
  t.charm       -1          # pin the tagged (charm -1) D_s- side
end

# Signal side: the recoil D_s+ (net charge +1)
alg_modeI.signal_side do |s|
  s.charged(pip: 3, pim: 2)  # 2 K_S0 -> 2 pi+ 2 pi-, plus the prompt pi+
  s.photons 3                # pi0 -> gamma gamma + the D_s*+ transition photon
  s.require_charge(1)        # net charge +1
  s.min_photon_angle 10.0    # photon opening angle > 10 deg to nearest track
  s.min_photon_energy 0.025  # barrel energy threshold 25 MeV
end

# Final 4C kinematic fit
alg_modeI.fit do |f|
  f.constrain_four_momentum  # 4C energy-momentum constraint
  f.chi2_cut 200             # chi2 < 200
end

# BOSS-side procedures that cannot be expressed in the tag DSL
alg_modeI
  .note(:tag_side_window, "per-mode M_BC and M_tag windows applied to the 16 hadronic D_s- tag modes; values are stored unconditionally and windowed in ROOT")
  .note(:background_veto, "for the D_s- -> pi- pi+ pi- tag mode only, pi+pi- pairs with invariant mass in [0.468, 0.528] GeV are vetoed (K_S0 window)")
  .note(:ks_reconstruction, "K_S0 built from pi+pi- with |M(pi+pi-) - m_K_S0| in [0.486, 0.510] GeV and flight distance > 2 sigma; among multiple candidates the longest-lived one(s) are retained")
  .note(:pi0_reconstruction, "pi0 built from gamma gamma with M in [0.115, 0.150] GeV and required to pass a 1C mass-constrained fit")
  .note(:pion_veto, "pions not originating from K_S0 / eta / eta' are rejected if p < 0.1 GeV/c")
  .note(:transition_candidate, "the D_s*+ transition photon (or pi0) is selected as the unused candidate minimising |DeltaE|")
  .note(:photon_selection, "photon barrel/endcap energy thresholds 25/50 MeV and TDC window [0,700] ns")

alg_modeI.apply   # no Selection argument for a tag analysis
alg_modeI.execute_on([ds_data, ds_incMC, exMC_modeI, exMC_modeII])

# ---------- Mode II: D_s+ -> K_S0 K+ pi0 pi0 ----------
alg_modeII = TagAnalysis.new("DsTagModeII")
alg_modeII.set_header(["DsTagModeIIAlg/DsTagModeII.h"])
          .set_constant({"ECMS" => [:double, 4.178]})
          .with_decay_card(decay_card_modeII)

# Tag side: same 16 hadronic D_s- modes
alg_modeII.tag_side(:Ds) do |t|
  t.mode_group(:hadronic)
  t.charm       -1
end

# Signal side: the recoil D_s+ (net charge +1)
alg_modeII.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)  # K+ and K_S0 -> pi+ pi-
  s.photons 5                        # 2 pi0 -> 4 gamma + the D_s*+ transition photon
  s.require_charge(1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Final 4C kinematic fit
alg_modeII.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_modeII
  .note(:tag_side_window, "per-mode M_BC and M_tag windows applied to the 16 hadronic D_s- tag modes; values are stored unconditionally and windowed in ROOT")
  .note(:background_veto, "for the D_s- -> pi- pi+ pi- tag mode only, pi+pi- pairs with invariant mass in [0.468, 0.528] GeV are vetoed (K_S0 window)")
  .note(:ks_reconstruction, "K_S0 built from pi+pi- with |M(pi+pi-) - m_K_S0| in [0.486, 0.510] GeV and flight distance > 2 sigma; among multiple candidates the longest-lived one(s) are retained")
  .note(:pi0_reconstruction, "pi0 built from gamma gamma with M in [0.115, 0.150] GeV passing a 1C fit; keeping the two candidates with the smallest chi2_1C")
  .note(:pion_veto, "pions not originating from K_S0 / eta / eta' are rejected if p < 0.1 GeV/c")
  .note(:transition_candidate, "the D_s*+ transition photon (or pi0) is selected as the unused candidate minimising |DeltaE|")
  .note(:photon_selection, "photon barrel/endcap energy thresholds 25/50 MeV and TDC window [0,700] ns")

alg_modeII.apply
alg_modeII.execute_on([ds_data, ds_incMC, exMC_modeI, exMC_modeII])