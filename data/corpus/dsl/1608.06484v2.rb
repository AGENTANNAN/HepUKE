### Dataset description ###
# 482 pb^-1 collected at sqrt(s) = 4.009 GeV with the BESIII detector.
data_4009    = DatasetManager.real_data.find("703_4009")
incMC_4009   = DatasetManager.inclusive_mc.find("703_4009")

### Decay cards for exclusive signal MC ###
decay_card_eta_enu = <<~DECAYCARD
    Decay psi(4040)
    1.0000  D_s+  D_s-                                     PHSP;
    Enddecay

    Decay D_s+
    1.0000  eta  e+  nu_e                                  ISGW2;
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

decay_card_etap_etapipi_enu = <<~DECAYCARD
    Decay psi(4040)
    1.0000  D_s+  D_s-                                     PHSP;
    Enddecay

    Decay D_s+
    1.0000  eta'  e+  nu_e                                 ISGW2;
    Enddecay

    Decay eta'
    1.0000  eta  pi+  pi-                                  PHSP;
    Enddecay

    Decay D_s-
    0.1000  K+  K-  pi-                                    PHSP;
    0.1000  phi  rho-                                      SVV_HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
    0.1000  K_S0  K+  pi-  pi-                             PHSP;
    0.1000  K_S0  K-  pi+  pi-                             PHSP;
    0.1000  K_S0  K-                                       PHSP;
    0.1000  pi+  pi-  pi-                                  PHSP;
    0.1000  eta  pi-                                       PHSP;
    0.1000  eta_tag_ppp  pi-                               PHSP;
    0.1000  eta_tag_grho  pi-                              PHSP;
    0.1000  eta  rho-                                      SVS;
    Enddecay

    Alias eta_tag_ppp eta'
    Alias eta_tag_grho eta'

    Decay eta_tag_ppp
    1.0000  eta  pi+  pi-                                  PHSP;
    Enddecay

    Decay eta_tag_grho
    1.0000  gamma  rho0                                    SVP_HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay phi
    1.0000  K+  K-                                         VSS;
    Enddecay

    Decay rho-
    1.0000  pi0  pi-                                       VSS;
    Enddecay

    Decay rho0
    1.0000  pi+  pi-                                       VSS;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                                   PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                                   PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                       PHSP;
    Enddecay

    End
DECAYCARD

decay_card_etap_grho_enu = <<~DECAYCARD
    Decay psi(4040)
    1.0000  D_s+  D_s-                                     PHSP;
    Enddecay

    Decay D_s+
    1.0000  eta'  e+  nu_e                                 ISGW2;
    Enddecay

    Decay eta'
    1.0000  gamma  rho0                                    SVP_HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay rho0
    1.0000  pi+  pi-                                       VSS;
    Enddecay

    Decay D_s-
    0.1000  K+  K-  pi-                                    PHSP;
    0.1000  phi  rho-                                      SVV_HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
    0.1000  K_S0  K+  pi-  pi-                             PHSP;
    0.1000  K_S0  K-  pi+  pi-                             PHSP;
    0.1000  K_S0  K-                                       PHSP;
    0.1000  pi+  pi-  pi-                                  PHSP;
    0.1000  eta  pi-                                       PHSP;
    0.1000  eta_tag_ppp  pi-                               PHSP;
    0.1000  eta_tag_grho  pi-                              PHSP;
    0.1000  eta  rho-                                      SVS;
    Enddecay

    Alias eta_tag_ppp eta'
    Alias eta_tag_grho eta'

    Decay eta_tag_ppp
    1.0000  eta  pi+  pi-                                  PHSP;
    Enddecay

    Decay eta_tag_grho
    1.0000  gamma  rho0                                    SVP_HELAMP 1.0 0.0 1.0 0.0;
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

    Decay pi0
    1.0000  gamma  gamma                                   PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                       PHSP;
    Enddecay

    End
DECAYCARD

exMC_eta_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Ds_eta_enu"
  config.related_dataset = data_4009
  config.events         = 500000
  config.decay_card     = decay_card_eta_enu
  config.cross_section  = :default
end
exMC_eta_enu.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_etap_etapipi_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Ds_etap_etapipi_enu"
  config.related_dataset = data_4009
  config.events         = 500000
  config.decay_card     = decay_card_etap_etapipi_enu
  config.cross_section  = :default
end
exMC_etap_etapipi_enu.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_etap_grho_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Ds_etap_grho_enu"
  config.related_dataset = data_4009
  config.events         = 500000
  config.decay_card     = decay_card_etap_grho_enu
  config.cross_section  = :default
end
exMC_etap_grho_enu.save_to_config(format: :yaml, file_path: 'temp_for_test')

datasets = [data_4009, incMC_4009]

ds_tag_modes = [
  :DstoKKPi,
  :DstoPhiRho,
  :DstoKsKPiPi,
  :DstoKsK,
  :DstoPiPiPi,
  :DstoPiEta,
  :DstoEtaPPi,
  :DstoEtaPPiGRho0,
  :DstoEtaRho
]

########################################################################
# Algorithm I : Ds+ -> eta e+ nu_e
########################################################################
alg_eta = TagAnalysis.new("DsTagEtaENu")
alg_eta.set_header(["DsTagEtaENuAlg/DsTagEtaENu.h"])
       .set_constant({ "ECMS" => [:double, 4.009] })
       .note(:tag_side_selection,
             "Ds- tag: 10 hadronic modes. Charged track: |cos(theta)|<0.93, " \
             "|Vz|<10 cm, |Vxy|<1 cm; PID pion (CL_pi>0, CL_pi>CL_K) or " \
             "kaon (CL_K>0, CL_K>CL_pi). Photon: E>25 MeV barrel |cos(theta)|<0.80 " \
             "or E>50 MeV endcap 0.86<|cos(theta)|<0.92, TDC in [0,700] ns, " \
             "angle to charged track > 10 deg. Intermediate windows: " \
             "0.115<M(gg)<0.150 (pi0), 0.510<M(gg)<0.570 (eta), " \
             "1.005<M(KK)<1.040 (phi), 0.570<M(pi0 pi-)<0.970 (rho-), " \
             "0.943<M(eta pi+ pi-)<0.973 and 0.932<M(gamma rho0)<0.980 " \
             "with 0.570<M(pi+pi-)<0.970 (eta'). K_S0: pipi pair Vz within " \
             "+-20 cm, 0.487<M(pi+pi-)<0.511, secondary vertex fit with " \
             "positive decay length. pi0/eta mass-constrained (1C).")
       .note(:tag_delta_e_window, "Delta E = E_ST - E_beam kept within +-3 sigma; " \
             "best tag candidate = minimum |Delta E|. Signal region: mBC within " \
             "(-4 sigma, 5 sigma) of the fitted peak. Applied downstream in ROOT.")
       .note(:electron_pid,
             "Signal electron: CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8. " \
             "Charge opposite to the tag Ds-.")
       .note(:extra_photon_energy_veto,
             "E_extragamma_max < 300 MeV: energy of the most energetic photon " \
             "not used in DT reconstruction must be below 0.300 GeV.")
       .note(:no_extra_track_veto,
             "No extra charged track beyond the DT candidates.")
       .note(:umiss_signal_region,
             "U_miss = E_miss - |p_miss| stored; signal region (-0.10, +0.12) GeV " \
             "applied at ROOT stage.")

alg_eta.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm -1
end

alg_eta.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_eta.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_eta.with_decay_card(decay_card_eta_enu).apply
alg_eta.execute_on(datasets + [exMC_eta_enu])

########################################################################
# Algorithm II : Ds+ -> eta' e+ nu_e, eta' -> eta(gg) pi+ pi-
########################################################################
alg_etap_ppp = TagAnalysis.new("DsTagEtaPEtaPiPiENu")
alg_etap_ppp.set_header(["DsTagEtaPEtaPiPiENuAlg/DsTagEtaPEtaPiPiENu.h"])
            .set_constant({ "ECMS" => [:double, 4.009] })
            .note(:tag_side_selection,
                  "Same ten Ds- tag modes and same tag-side track/photon/" \
                  "intermediate-state windows as the eta e+ nu_e algorithm.")
            .note(:electron_pid,
                  "CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8; charge opposite " \
                  "to the tag Ds-.")
            .note(:etap_mass_window,
                  "0.943 < M(eta pi+ pi-) < 0.973 GeV/c^2; eta -> gamma gamma with " \
                  "0.510 < M(gg) < 0.570 GeV/c^2 and 1C mass constraint applied.")
            .note(:extra_photon_energy_veto, "E_extragamma_max < 300 MeV.")
            .note(:no_extra_track_or_pi0_veto,
                  "No extra charged track and no extra pi0 beyond the DT candidates.")
            .note(:umiss_signal_region,
                  "U_miss signal region (-0.10, +0.12) GeV applied in ROOT.")

alg_etap_ppp.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm -1
end

alg_etap_ppp.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(ep: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_etap_ppp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_etap_ppp.with_decay_card(decay_card_etap_etapipi_enu).apply
alg_etap_ppp.execute_on(datasets + [exMC_etap_etapipi_enu])

########################################################################
# Algorithm III : Ds+ -> eta' e+ nu_e, eta' -> gamma rho0
########################################################################
alg_etap_grho = TagAnalysis.new("DsTagEtaPGRhoENu")
alg_etap_grho.set_header(["DsTagEtaPGRhoENuAlg/DsTagEtaPGRhoENu.h"])
             .set_constant({ "ECMS" => [:double, 4.009] })
             .note(:tag_side_selection,
                   "Same ten Ds- tag modes and same tag-side selection as the " \
                   "other two algorithms.")
             .note(:electron_pid,
                   "CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8; charge opposite " \
                   "to the tag Ds-.")
             .note(:etap_mass_window,
                   "0.932 < M(gamma rho0) < 0.980 GeV/c^2 with " \
                   "0.570 < M(pi+ pi-) < 0.970 GeV/c^2 for the rho0.")
             .note(:extra_photon_energy_veto, "E_extragamma_max < 300 MeV.")
             .note(:no_extra_track_or_pi0_veto,
                   "No extra charged track and no extra pi0 beyond the DT candidates.")
             .note(:umiss_signal_region,
                   "U_miss signal region (-0.08, +0.10) GeV applied in ROOT.")

alg_etap_grho.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm -1
end

alg_etap_grho.signal_side do |s|
  s.photons 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(ep: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_etap_grho.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_etap_grho.with_decay_card(decay_card_etap_grho_enu).apply
alg_etap_grho.execute_on(datasets + [exMC_etap_grho_enu])
