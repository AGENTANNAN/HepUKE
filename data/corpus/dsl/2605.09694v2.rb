# ============================================================
# Dataset preparation
# ============================================================
# Ds Ds* data collected at CMS energies 4.128 - 4.226 GeV (7.33 fb^-1).
# Use the ~4.178 GeV working point as the anchor dataset.
data_4178 = DatasetManager.real_data.find("705_4180")
incMC_4178 = DatasetManager.inclusive_mc.find("705_4180")

# --------------------------------------------------------------------
# Decay card for signal mode I: Ds+ -> K_S0 K_S0 pi+ pi0
# --------------------------------------------------------------------
decay_card_modeI = <<~DECAYCARD
  Decay D_s+
  1.0000  K_S0  K_S0  pi+  pi0     PHSP;
  Enddecay
  CDecay D_s-

  Decay K_S0
  1.0000  pi+  pi-                 PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma gamma              PHSP;
  Enddecay

  End
DECAYCARD

# --------------------------------------------------------------------
# Decay card for signal mode II: Ds+ -> K_S0 K+ pi0 pi0
# --------------------------------------------------------------------
decay_card_modeII = <<~DECAYCARD
  Decay D_s+
  1.0000  K_S0  K+  pi0  pi0       PHSP;
  Enddecay
  CDecay D_s-

  Decay K_S0
  1.0000  pi+  pi-                 PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma gamma              PHSP;
  Enddecay

  End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Ds_to_KsKsPiPi0"
  config.related_dataset = data_4178
  config.events         = 500000
  config.decay_card     = decay_card_modeI
  config.cross_section  = :default
end
exMC_modeI.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Ds_to_KsKPi0Pi0"
  config.related_dataset = data_4178
  config.events         = 500000
  config.decay_card     = decay_card_modeII
  config.cross_section  = :default
end
exMC_modeII.save_to_config(format: :yaml, file_path: 'temp_for_test')

# ============================================================
# Event selection: DTag double-tag reconstruction
# ST: Ds- reconstructed in 16 hadronic modes.
# DT: full hadronic reconstruction of the signal Ds+ decay,
#     including the transition gamma or pi0 from Ds*+.
# ============================================================

# --------------------------------------------------------
# Mode I: Ds+ -> K_S0 K_S0 pi+ pi0
# --------------------------------------------------------
alg_modeI = TagAnalysis.new("DsToKsKsPiPi0")
alg_modeI.set_header(["DsToKsKsPiPi0Alg/DsToKsKsPiPi0.h"])
         .set_constant({ "ECMS" => [:double, 4.178] })
         .with_decay_card(decay_card_modeI)

# Single-tag anti-Ds- reconstructed in the sixteen hadronic modes
alg_modeI.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKKPiPi0,
          :DstoPiPiPi,
          :DstoKsK,
          :DstoKsKPi0,
          :DstoKPiPi,
          :DstoKsKsPi,
          :DstoKsKPiPi,
          :DstoKsKPiPi,
          :DstoPiEta,
          :DstoPiPi0Eta,
          :DstoPiEtaPr,
          :DstoPiEtaPr,
          :DstoRhoEta,
          :DstoRhoEta,
          :DstoPiPiPiEta
  t.charm -1
end

# Signal side: 2 K_S0 (-> pi+pi-) + pi+ + pi0 (-> gg) + transition gamma
# Charged tracks: 3 pi+ and 2 pi- (net charge +1)
# Photons: 2 for pi0 + 1 transition gamma = >=3
alg_modeI.signal_side do |s|
  s.charged(pip: 3, pim: 2)
  s.require_charge 1
  s.photons 3
  s.min_photon_angle 10.0
end

alg_modeI.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_modeI.note(:ks_reconstruction,
               "Two K_S0 candidates reconstructed from pi+ pi- pairs with a secondary vertex fit; " \
               "|M(pi+pi-) - m_Ks| in (0.486, 0.510) GeV; flight distance > 2 sigma. " \
               "If more than two Ks candidates, keep the two with the longest decay lengths.")
         .note(:pi0_reconstruction,
               "pi0 reconstructed from gamma-gamma with |M_gg - m_pi0| in (0.115, 0.150) GeV and a 1C kinematic fit. " \
               "Photon requirements: E_barrel > 25 MeV, E_endcap > 50 MeV, TDC in [0,700] ns, opening angle to nearest track > 10 deg.")
         .note(:soft_pion_veto,
               "Momentum of any pi not from K_S0/eta/eta' must be > 0.1 GeV/c to veto soft pions from D*+.")
         .note(:transition_photon_selection,
               "Transition gamma or pi0 from Ds*+ chosen as the unused candidate minimizing " \
               "|Delta E| = |Ecm - E_tag - E_recoil(gamma/pi0,Ds) - E_gamma(pi0)|.")
         .note(:tag_side,
               "ST Ds- reconstructed in 16 hadronic decay modes; per-tag M_BC and M_tag windows " \
               "(see Table 1 and Table 2 of the paper) applied in ROOT stage.")
         .note(:peaking_bkg_veto_3pi_tag,
               "For the tag Ds- -> pi- pi+ pi-, veto any pi+pi- combination with mass in (0.468, 0.528) GeV " \
               "to suppress D- -> K_S0 pi- peaking background.")

alg_modeI.apply
alg_modeI.execute_on([data_4178, incMC_4178, exMC_modeI])

# --------------------------------------------------------
# Mode II: Ds+ -> K_S0 K+ pi0 pi0
# --------------------------------------------------------
alg_modeII = TagAnalysis.new("DsToKsKPi0Pi0")
alg_modeII.set_header(["DsToKsKPi0Pi0Alg/DsToKsKPi0Pi0.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })
          .with_decay_card(decay_card_modeII)

alg_modeII.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKKPiPi0,
          :DstoPiPiPi,
          :DstoKsK,
          :DstoKsKPi0,
          :DstoKPiPi,
          :DstoKsKsPi,
          :DstoKsKPiPi,
          :DstoKsKPiPi,
          :DstoPiEta,
          :DstoPiPi0Eta,
          :DstoPiEtaPr,
          :DstoPiEtaPr,
          :DstoRhoEta,
          :DstoRhoEta,
          :DstoPiPiPiEta
  t.charm -1
end

# Signal side: K_S0(->pi+pi-) + K+ + 2 pi0 (each -> gamma gamma) + transition gamma
# Charged tracks: 1 K+, 1 pi+, 1 pi- (net charge +1)
# Photons: 4 (from 2 pi0) + 1 transition gamma = >=5
alg_modeII.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.photons 5
  s.min_photon_angle 10.0
end

alg_modeII.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_modeII.note(:ks_reconstruction,
                "K_S0 from pi+pi- with secondary vertex fit; |M(pi+pi-) - m_Ks| in (0.486, 0.510) GeV; flight distance > 2 sigma; " \
                "if more than one Ks candidate, keep the one with the longest decay length.")
          .note(:pi0_reconstruction,
                "Two pi0's reconstructed from gamma-gamma pairs with mass window (0.115, 0.150) GeV and 1C kinematic fit; " \
                "if more than two pi0 candidates, keep the two with smallest chi^2_1C.")
          .note(:soft_pion_veto,
                "Momentum of any pi not from K_S0/eta/eta' must be > 0.1 GeV/c.")
          .note(:transition_photon_selection,
                "Transition gamma or pi0 from Ds*+ chosen as the unused candidate minimizing |Delta E|.")
          .note(:tag_side,
                "ST Ds- reconstructed in 16 hadronic decay modes; per-tag M_BC and M_tag windows applied in ROOT.")
          .note(:peaking_bkg_veto_3pi_tag,
                "For tag Ds- -> pi- pi+ pi-, veto any pi+pi- with mass in (0.468, 0.528) GeV to remove K_S0 pi- contamination.")

alg_modeII.apply
alg_modeII.execute_on([data_4178, incMC_4178, exMC_modeII])
