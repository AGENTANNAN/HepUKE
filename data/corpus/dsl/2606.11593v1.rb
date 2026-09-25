### Dataset description ###
# Energy-scan samples 4.128 - 4.226 GeV (total 7.33 fb^-1), used for Ds*+ Ds- production
data_points = [
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230")
]

inc_mc_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

########################################################################
# Common Ds- tag mode list (eleven hadronic modes)
########################################################################
DS_TAG_MODES = [
  :DstoKsK, :DstoKKPi, :DstoKKPiPi0,
  :DstoKsKPiPi, :DstoKsKPiPi,           # KsK-pi+pi- and KsK+pi-pi-
  :DstoPiPiPi,
  :DstoPiEta, :DstoPiPi0Eta,
  :DstoPiEPPiPiEta, :DstoPiEtaPrimeGammaRho0,
  :DstoKPiPi
].uniq

########################################################################
# Signal decay card: e+e- -> Ds*+ Ds- ; Ds*+ -> gamma/pi0 Ds+ ; Ds+ -> tau+ nu ; tau -> ...
# One card covers all four tau decays because we use PHSP everywhere; the four
# Algorithm objects share it.
########################################################################
decay_card_signal = <<~DECAYCARD
    Decay vpho
    0.5000  D_s*+ D_s-                                PHSP;
    0.5000  D_s*- D_s+                                PHSP;
    Enddecay

    Decay D_s*+
    0.9350  gamma D_s+                                VSP_PWAVE;
    0.0650  pi0 D_s+                                  PHSP;
    Enddecay

    Decay D_s*-
    0.9350  gamma D_s-                                VSP_PWAVE;
    0.0650  pi0 D_s-                                  PHSP;
    Enddecay

    Decay D_s+
    1.0000  tau+ nu_tau                               SLN;
    Enddecay

    Decay D_s-
    0.0500  K_S0 K-                                   PHSP;
    Enddecay

    Decay tau+
    0.2500  e+     nu_e     anti-nu_tau               PHSP;
    0.2500  mu+    nu_mu    anti-nu_tau               PHSP;
    0.2500  pi+             anti-nu_tau               PHSP;
    0.2500  pi+ pi0         anti-nu_tau               PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+ pi-                                   PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                               PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "DsToTauNu_signal"
  config.events        = 500000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

########################################################################
# Helper: common tag-side declaration + notes shared by all four tau modes
########################################################################
add_common_notes = lambda do |a|
  a.note(:tag_side_selection,
         "eleven Ds- hadronic tag modes: K_S0 K-, K+K- pi-, K+K- pi- pi0, K_S0 K- pi+ pi-, " \
         "K_S0 K+ pi- pi-, pi+ pi- pi-, pi- eta, pi- pi0 eta, pi- eta'(pi+pi-eta), " \
         "pi- eta'(gamma rho0), K- pi+ pi-. Charged tracks (non-KS) require |Vz| < 10 cm, " \
         "|Vxy| < 1 cm and |cos(theta)| < 0.93. K/pi PID via L(K) vs L(pi) using MDC+TOF. " \
         "KS -> pi+pi- via vertex fit chi2 < 100, M(pi+pi-) in [0.487,0.511] GeV/c^2 and " \
         "decay length > 2*sigma. Photons: E > 25 (50) MeV in barrel (endcap); " \
         "|cos(theta)| in barrel < 0.80 or endcap 0.86<|c|<0.92; angle-to-track > 10 deg; " \
         "EMC time in [0,700] ns. pi0/eta from gamma-gamma with |M-m|<mass windows and 1C fit. " \
         "rho0 in [0.500,0.970] GeV/c^2; eta' -> pi+pi-eta in [0.946,0.970], eta' -> gamma rho0 in " \
         "[0.940,0.976] GeV/c^2. Momentum > 0.1 GeV/c required for pion and gamma from eta'-> gamma rho0. " \
         "The Ds- -> K_S0 K- and K- pi+ pi- overlap is removed by vetoing M(pi+pi-) in [0.480,0.515] GeV/c^2 " \
         "for the K- pi+ pi- mode.")
    .note(:ds_recoil_mass,
         "Recoil mass Mrec against Ds- required in [2.050, 2.195] GeV/c^2. If multiple ST Ds- " \
         "candidates, the one with the minimum |Mrec - m_Ds*+| is retained for each mode.")
    .note(:ds_mass_window,
         "Reconstructed Ds- invariant mass required within the mode-dependent window " \
         "(approximately +/-3 sigma around the Ds- mass, see Table I of the paper).")
    .note(:transition_photon_pi0,
         "In the presence of the best ST Ds-, a transition gamma or pi0 from Ds*+ decay is " \
         "reconstructed with the same criteria as tag-side gamma/pi0. Delta_E = Ecm - E_ST - " \
         "E_miss' - E_gamma(pi0) must be within [-0.2, 0.2] GeV. Multiple gamma/pi0 candidates: " \
         "the one with minimum |Delta_E| is kept.")
end

########################################################################
# Mode 1: tau+ -> e+ nu_e anti-nu_tau  (ST tag + missing Ds*+ side)
########################################################################
alg_e = TagAnalysis.new("DsToTauNu_TauE")
alg_e.set_header(["DsToTauNu_TauEAlg/DsToTauNu_TauE.h"])
     .set_constant({"ECMS" => [:double, 4.178]})
     .with_decay_card(decay_card_signal)

alg_e.tag_side(:Ds) { |t| t.modes(*DS_TAG_MODES); t.charm(-1) }

alg_e.signal_side do |s|
  s.photons 1..48                        # transition gamma from Ds*+ (also allow pi0-mode via extra photons)
  s.charged(ep: 1)                       # e+ from tau+ decay
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :X0, mass: nil               # missing neutrino system (mass ~0)
end

alg_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.store_fitted_momenta
end

add_common_notes.call(alg_e)
alg_e.note(:tau_e_selection,
           "Additional good track (opposite charge to ST Ds-) identified as positron: L(e) > 0.001, " \
           "L(e)/(L(e)+L(pi)+L(K)) > 0.8; p(e+) > 0.2 GeV/c; E_EMC/p(e+) > 0.8.")

alg_e.apply
alg_e.execute_on(data_points + inc_mc_points + exMC_signal)

########################################################################
# Mode 2: tau+ -> mu+ nu_mu anti-nu_tau
########################################################################
alg_mu = TagAnalysis.new("DsToTauNu_TauMu")
alg_mu.set_header(["DsToTauNu_TauMuAlg/DsToTauNu_TauMu.h"])
      .set_constant({"ECMS" => [:double, 4.178]})
      .with_decay_card(decay_card_signal)

alg_mu.tag_side(:Ds) { |t| t.modes(*DS_TAG_MODES); t.charm(-1) }

alg_mu.signal_side do |s|
  s.photons 1..48
  s.charged(mup: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :X0, mass: nil
end

alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.store_fitted_momenta
end

add_common_notes.call(alg_mu)
alg_mu.note(:tau_mu_selection,
            "Additional good track (opposite charge to ST Ds-) identified as muon: L(mu) > 0.001, " \
            "L(mu) > L(K), L(mu) > L(e); E_EMC in [0.1, 0.3] GeV; p(mu+) > 0.5 GeV/c; MUC hit " \
            "depth requirements are momentum- and polar-angle dependent (Table II of paper).")

alg_mu.apply
alg_mu.execute_on(data_points + inc_mc_points + exMC_signal)

########################################################################
# Mode 3: tau+ -> pi+ anti-nu_tau
########################################################################
alg_pi = TagAnalysis.new("DsToTauNu_TauPi")
alg_pi.set_header(["DsToTauNu_TauPiAlg/DsToTauNu_TauPi.h"])
      .set_constant({"ECMS" => [:double, 4.178]})
      .with_decay_card(decay_card_signal)

alg_pi.tag_side(:Ds) { |t| t.modes(*DS_TAG_MODES); t.charm(-1) }

alg_pi.signal_side do |s|
  s.photons 1..48
  s.charged(pip: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :X0, mass: nil
end

alg_pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.store_fitted_momenta
end

add_common_notes.call(alg_pi)
alg_pi.note(:tau_pi_selection,
            "Additional good track (opposite charge to ST Ds-) identified as pi+: L(pi) > L(K); " \
            "no pi0 candidate allowed (to suppress tau_rho contamination); E_EMC/p(pi+) < 0.9; " \
            "total energy of residual good photons < 0.3 GeV; cos(theta_miss) < 0.9.")

alg_pi.apply
alg_pi.execute_on(data_points + inc_mc_points + exMC_signal)

########################################################################
# Mode 4: tau+ -> pi+ pi0 anti-nu_tau  (with pi+pi0 in rho+ region)
########################################################################
alg_rho = TagAnalysis.new("DsToTauNu_TauRho")
alg_rho.set_header(["DsToTauNu_TauRhoAlg/DsToTauNu_TauRho.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
       .with_decay_card(decay_card_signal)

alg_rho.tag_side(:Ds) { |t| t.modes(*DS_TAG_MODES); t.charm(-1) }

alg_rho.signal_side do |s|
  s.photons 3..48                        # transition gamma from Ds*+ + pi0 from tau+ (>=2 photons for pi0)
  s.charged(pip: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :X0, mass: nil
end

alg_rho.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # constrain the tau-side pi0
  f.chi2_cut 200
  f.store_fitted_momenta
end

add_common_notes.call(alg_rho)
alg_rho.note(:tau_rho_selection,
             "Additional good track identified as pi+ (as tau_pi); the tau-side pi0 is reconstructed " \
             "from a good photon pair, best pi0 by minimum chi2 of the 1C kinematic fit; total energy " \
             "of additional good photons < 0.1 GeV. |M(pi+ pi0) - m_rho+| < 0.2 GeV/c^2 to enforce " \
             "the rho+ region. Missing mass squared M'^2_miss required in [3.82, 3.98] GeV^2/c^4.")

alg_rho.apply
alg_rho.execute_on(data_points + inc_mc_points + exMC_signal)
