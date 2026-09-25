# BESIII: D+ -> K- pi+ e+ nu_e (D_e4 semileptonic decay)
# Data: 2.93 fb^-1 at sqrt(s) = 3.773 GeV (psi(3770)); tagged-D technique.
# Tag side: six hadronic D- modes. Signal side: K- pi+ e+ + missing nu_e.

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process D+ -> K- pi+ e+ nu_e
# (generated with the decay intensity distribution determined by PWA;
#  the D- opposite the signal decays inclusively as in generic MC)
decay_card_D_e4 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                              PHSP;
  Enddecay

  Decay D+
  1.0000 K- pi+ e+ nu_e                     PHOTOS PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                         PHSP;
  Enddecay

  End
DECAYCARD

exMC_D_e4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_to_Kmpip_ep_nue"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_D_e4
  config.cross_section   = :default
end
exMC_D_e4.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — TagAnalysis (single tag + missing neutrino) ###

alg = TagAnalysis.new("DpToKmpipEpNue")
alg.set_header(["DpToKmpipEpNueAlg/DpToKmpipEpNue.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })     # sqrt(s) = 3.773 GeV
   .with_decay_card(decay_card_D_e4)
   .note(:tag_side_selection,
         "D- tag reconstructed in six hadronic channels with the largest " \
         "branching fractions. Good track: |Vz| < 10 cm, |Vr| < 1 cm, " \
         "|cos(theta)| < 0.93. Photon (for pi0): E > 25 MeV barrel " \
         "(|cos(theta)| < 0.80) or E > 50 MeV endcap (0.86 < |cos(theta)| < 0.92), " \
         "shower time within 700 ns after the beam collision. " \
         "pi/K PID from combined dE/dx and TOF likelihoods: P(K) > P(pi) for kaons, " \
         "P(pi) > P(K) for pions. pi0: 0.115 < M(gamma gamma) < 0.150 GeV/c^2 with a " \
         "1-C mass-constrained fit chi2 < 200; events with both photons in the " \
         "endcap rejected. K_S0: secondary vertex from oppositely charged tracks " \
         "(no good-track or PID requirement on the daughters), " \
         "0.487 < M(pi+ pi-) < 0.511 GeV/c^2, closest approach within 20 cm of the " \
         "IP along the beam direction and |cos(theta)| < 0.93.")
   .note(:tag_delta_e_window,
         "Tag Delta E = E_D - E_beam required consistent with zero within " \
         "approximately twice the experimental resolution (mode-dependent, " \
         "applied downstream in ROOT). At most one candidate per tag mode per " \
         "charge is accepted; for multiple candidates the one with the smallest " \
         "|Delta E| is chosen.")
   .note(:tag_mbc_sideband,
         "The tagged decay yields are extracted from a fit to the M_BC distribution " \
         "(signal shape from reconstructed MC, background shape from the ARGUS " \
         "function); the M_BC sideband of the data provides the combinatorial " \
         "background estimate.")
   .note(:signal_side_track_multiplicity,
         "Exactly three signal-side tracks satisfying the good-track criteria, " \
         "identified as K-, pi+ and e+. Electron PID: combined dE/dx, TOF and EMC " \
         "likelihoods with P2(e)/(P2(K)+P2(pi)+P2(e)) > 0.8 and P2(e) > 0.001; " \
         "EMC energy of the electron candidate > 80% of the MDC track momentum.")
   .note(:missing_neutrino,
         "E_miss and p_miss of the neutrino are reconstructed from energy and " \
         "momentum conservation. Undetected massive-particle background suppressed " \
         "by |U_miss| < 0.04 GeV with U_miss = E_miss - |p_miss|; neutrino-less " \
         "decays suppressed by E_miss > 0.04 GeV.")
   .note(:extra_shower_veto,
         "Background from events containing pi0 suppressed by requiring that no " \
         "unassociated EMC shower has an energy deposition above 0.25 GeV " \
         "(clusters separated by more than 15 degrees from the closest charged " \
         "tracks are considered).")
   .note(:d0_crossfeed_veto,
         "Cross-feed from e+e- -> D0 D0bar rejected for the tag modes " \
         "K_S0 pi- pi- pi+, K_S0 pi- pi0 and K+ pi- pi- pi0: a purely hadronic D0 " \
         "decay is reconstructed from all event tracks and the event is rejected if " \
         "1.860 < M_BC(D0) < 1.875 GeV/c^2 and |Delta E(D0)| < 0.01 GeV.")
   .note(:pwa_component_fractions,
         "Partial wave analysis (unbinned maximum-likelihood fit in the five " \
         "kinematic variables m^2, q^2, cos theta_K, cos theta_e, chi) is performed " \
         "downstream; the K*(892)^0-dominant solution plus a K pi S-wave component " \
         "(f_S = (6.05 +- 0.22 +- 0.18)%) is used to weight the PWA signal MC and to " \
         "obtain the selection efficiency eps_tag,sig.")
   .note(:helicity_form_factors,
         "Model-independent helicity-basis form-factor products extracted with the " \
         "projective weighting technique in the K*-dominated region " \
         "(0.8 < m_Kpi < 1.0 GeV/c^2; 10 q^2 bins over 0 < q^2 < 1.0 GeV^2/c^4 and " \
         "100 angular bins per q^2 bin) using PHSP signal MC.")

# --- Tag side: D- single tag in six hadronic modes ---
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,        # D- -> K+ pi- pi-
          :DptoKPiPiPi0,     # D- -> K+ pi- pi- pi0
          :DptoKsPi,         # D- -> K_S0 pi-
          :DptoKsPiPi0,      # D- -> K_S0 pi- pi0
          :DptoKsPiPiPi,     # D- -> K_S0 pi- pi- pi+
          :DptoKKPi          # D- -> K+ K- pi-
  t.charm -1                 # tagged D- meson
  t.window :mBC, min: 1.863, max: 1.877   # M_BC signal region
end

# --- Signal side: K- pi+ e+ from the tracks not used by the tag, plus missing nu_e ---
alg.signal_side do |s|
  s.charged(km: 1, pip: 1, ep: 1)          # exactly three signal-side tracks
  s.require_charge 1                       # -1 (K-) + 1 (pi+) + 1 (e+) = +1
  s.missing :nu_e                          # massless neutrino (semileptonic)
  s.min_photon_angle 10.0
end

# --- Kinematic fit: 4C total four-momentum (tag + K- pi+ e+ + nu = CMS) ---
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200                           # loose in BOSS; tight cut applied in ROOT
end

alg.apply
root_files = alg.execute_on([psi3770_data, psi3770_incMC, exMC_D_e4])
