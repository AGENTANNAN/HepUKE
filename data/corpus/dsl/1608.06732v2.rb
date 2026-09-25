### Dataset description ###
# 482 pb^-1 collected at sqrt(s) = 4.009 GeV with the BESIII detector.
data_4009    = DatasetManager.real_data.find("703_4009")
incMC_4009   = DatasetManager.inclusive_mc.find("703_4009")

### Decay cards for exclusive signal MC ###
# Ds+ -> mu+ nu_mu; the tag Ds- decays generically over the nine hadronic modes.
decay_card_mu_nu = <<~DECAYCARD
    Decay psi(4040)
    1.0000  D_s+  D_s-                                     PHSP;
    Enddecay

    Decay D_s+
    1.0000  mu+  nu_mu                                     PHSP;
    Enddecay

    Decay D_s-
    0.1000  K+  K-  pi-                                    PHSP;
    0.1000  phi  rho-                                      SVV_HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
    0.1000  K_S0  K+  pi-  pi-                             PHSP;
    0.1000  K_S0  K-  pi+  pi-                             PHSP;
    0.1000  K_S0  K-                                       PHSP;
    0.1000  pi+  pi-  pi-                                  PHSP;
    0.1000  eta  pi-                                       PHSP;
    0.1000  eta'  pi-                                      PHSP;
    0.1000  eta'  pi-                                      PHSP;
    0.1000  eta  rho-                                      SVS;
    Enddecay

    Decay phi
    1.0000  K+  K-                                         VSS;
    Enddecay

    Decay rho-
    1.0000  pi0  pi-                                       VSS;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                                   PHSP;
    Enddecay

    Decay eta'
    0.5000  eta  pi+  pi-                                  PHSP;
    0.5000  gamma  rho0                                    SVP_HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay rho0
    1.0000  pi+  pi-                                       VSS;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                                   PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                       PHSP;
    Enddecay

    End
DECAYCARD

# Ds+ -> tau+ nu_tau, tau+ -> pi+ nu_tau (tauonic mode).
decay_card_tau_nu = <<~DECAYCARD
    Decay psi(4040)
    1.0000  D_s+  D_s-                                     PHSP;
    Enddecay

    Decay D_s+
    1.0000  tau+  nu_tau                                   PHSP;
    Enddecay

    Decay tau+
    1.0000  pi+  nu_tau                                    PHSP;
    Enddecay

    Decay D_s-
    0.1000  K+  K-  pi-                                    PHSP;
    0.1000  phi  rho-                                      SVV_HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
    0.1000  K_S0  K+  pi-  pi-                             PHSP;
    0.1000  K_S0  K-  pi+  pi-                             PHSP;
    0.1000  K_S0  K-                                       PHSP;
    0.1000  pi+  pi-  pi-                                  PHSP;
    0.1000  eta  pi-                                       PHSP;
    0.1000  eta'  pi-                                      PHSP;
    0.1000  eta'  pi-                                      PHSP;
    0.1000  eta  rho-                                      SVS;
    Enddecay

    Decay phi
    1.0000  K+  K-                                         VSS;
    Enddecay

    Decay rho-
    1.0000  pi0  pi-                                       VSS;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                                   PHSP;
    Enddecay

    Decay eta'
    0.5000  eta  pi+  pi-                                  PHSP;
    0.5000  gamma  rho0                                    SVP_HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay rho0
    1.0000  pi+  pi-                                       VSS;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                                   PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                       PHSP;
    Enddecay

    End
DECAYCARD

exMC_mu_nu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Ds_mu_nu"
  config.related_dataset = data_4009
  config.events          = 500000
  config.decay_card      = decay_card_mu_nu
  config.cross_section   = :default
end

exMC_tau_nu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Ds_tau_nu"
  config.related_dataset = data_4009
  config.events          = 500000
  config.decay_card      = decay_card_tau_nu
  config.cross_section   = :default
end

datasets = [data_4009, incMC_4009]

# Nine hadronic Ds- single-tag modes: K_S0 K-, K+ K- pi-, K+ K- pi- pi0,
# K_S0 K+ pi- pi-, pi+ pi- pi-, eta pi- (eta->gg), pi0 eta pi-,
# eta' pi- (eta'->pi+ pi- eta) and eta' pi- (eta'->pi+ pi- gamma).
ds_tag_modes = [
  :DstoKsK,
  :DstoKKPi,
  :DstoKKPiPi0,
  :DstoKsKPiPi,
  :DstoPiPiPi,
  :DstoPiEta,
  :DstoPiPi0Eta,
  :DstoEtaPPiPiEta,
  :DstoEtaPRhoGam
]

########################################################################
# Algorithm I : Ds+ -> mu+ nu_mu
########################################################################
alg_mu = TagAnalysis.new("DsTagMuNu")
alg_mu.set_header(["DsTagMuNuAlg/DsTagMuNu.h"])
      .set_constant({ "ECMS" => [:double, 4.009] })
      .note(:tag_side_selection,
            "Ds- single tag from nine hadronic modes. Charged track: |cos(theta)|<0.93, " \
            "|Vz|<10 cm along the beam and Vr<1 cm perpendicular; PID requires " \
            "CL_pi>0 and CL_pi>CL_K for pions, CL_K>0 and CL_K>CL_pi for kaons. " \
            "Photon: E>25 MeV barrel (|cos(theta)|<0.80) or E>50 MeV endcap " \
            "(0.86<|cos(theta)|<0.92), TDC within [0,700] ns, angle to any charged " \
            "track >10 deg. Intermediate windows: 0.115<M(gg)<0.150 GeV/c^2 (pi0), " \
            "0.510<M(gg)<0.570 GeV/c^2 (eta), 0.943<M(pi+pi-eta)<0.973 GeV/c^2 and " \
            "0.932<M(gamma rho0)<0.980 GeV/c^2 with 0.570<M(pi+pi-)<0.970 GeV/c^2 " \
            "(eta'). K_S0: pi+pi- pair with |Vz|<20 cm and 0.487<M(pi+pi-)<0.511 " \
            "GeV/c^2, secondary-vertex fit with decay length at least twice the " \
            "vertex resolution. pi0/eta mass-constrained (1C).")
      .note(:tag_mbc_window,
            "Ds- tag candidates selected in 1.962 < M_BC < 1.982 GeV/c^2; the M_BC " \
            "distribution is fitted (signal shape from MC plus a background shape) " \
            "to extract the tag yields, totalling 15127 +- 321 tags.")
      .note(:signal_side_selection,
            "Signal side: exactly one good charged track with charge opposite to the " \
            "tag Ds-, treated as a muon with NO PID requirement; the most energetic " \
            "neutral shower not associated with the tag must carry less than 300 MeV.")
      .note(:missing_mass_signal_region,
            "MM^2 = (E_beam - E_mu+)^2/c^4 - (-p_Ds- - p_mu+)^2/c^2 is stored for every " \
            "event; the signal region -0.15 < MM^2 < 0.20 (GeV/c^2)^2 is applied " \
            "downstream in ROOT. Detection efficiency 91.4 +- 0.5%.")
      .note(:no_kinematic_constraint,
            "Because of the undetected neutrino no kinematic fit is imposed; only the " \
            "four-momentum conservation constraint is recorded for the fit block.")

alg_mu.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm -1
end

alg_mu.signal_side do |s|
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu.with_decay_card(decay_card_mu_nu).apply
alg_mu.execute_on(datasets + [exMC_mu_nu])

########################################################################
# Algorithm II : Ds+ -> tau+ nu_tau, tau+ -> pi+ nu_tau
########################################################################
alg_tau = TagAnalysis.new("DsTagTauNu")
alg_tau.set_header(["DsTagTauNuAlg/DsTagTauNu.h"])
       .set_constant({ "ECMS" => [:double, 4.009] })
       .note(:tag_side_selection,
             "Identical Ds- single-tag selection and nine tag modes as the " \
             "Ds+ -> mu+ nu_mu algorithm.")
       .note(:tag_mbc_window,
             "1.962 < M_BC < 1.982 GeV/c^2 for the tag Ds-, as in the muonic algorithm.")
       .note(:signal_side_selection,
             "Signal side: exactly one good charged track with charge opposite to the " \
             "tag Ds-, treated as a pion with NO PID requirement, plus the most " \
             "energetic neutral shower not associated with the tag below 300 MeV.")
       .note(:missing_mass_signal_region,
             "For Ds+ -> tau+ nu_tau with tau+ -> pi+ nu_tau there are two undetected " \
             "neutrinos; MM^2 = (E_beam - E_pi+)^2/c^4 - (-p_Ds- - p_pi+)^2/c^2 is " \
             "stored and the signal region -0.15 < MM^2 < 0.20 (GeV/c^2)^2 is applied " \
             "in ROOT. Detection efficiency 41.0 +- 0.3%.")

alg_tau.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm -1
end

alg_tau.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
  s.missing :nu_tau
end

alg_tau.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_tau.with_decay_card(decay_card_tau_nu).apply
alg_tau.execute_on(datasets + [exMC_tau_nu])
