# arXiv: 1901.02133v2
# Measurement of BFs for Ds+ -> eta(') e+ nu_e and form factors
# BESIII Collaboration
# Double-tag (ST+DT) method at sqrt(s) = 4.178 GeV, 3.19 fb^-1

# ============================================================================
# Dataset preparation
# ============================================================================

ds_4180 = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# ============================================================================
# TagAnalysis: Ds tag + semileptonic signal
# ST: Ds- reconstructed via 14 hadronic modes
# DT: Ds+ -> eta(') e+ nu_e
# ============================================================================

# --- Semileptonic decay cards for exclusive MC ---
decay_card_eta_enu = <<~DECAYCARD
  Decay psi(4260)
  1.000  D_s- D_s+               PHSP;
  Enddecay

  Decay D_s-
  1.000  K+ K- pi-               PHSP;
  Enddecay

  Decay D_s+
  1.000  eta e+ nu_e             SLN;
  Enddecay

  Decay eta
  1.000  gamma gamma             PHSP;
  Enddecay

  End
DECAYCARD

decay_card_etap_enu = <<~DECAYCARD
  Decay psi(4260)
  1.000  D_s- D_s+               PHSP;
  Enddecay

  Decay D_s-
  1.000  K+ K- pi-               PHSP;
  Enddecay

  Decay D_s+
  1.000  eta' e+ nu_e            SLN;
  Enddecay

  Decay eta'
  1.000  eta pi+ pi-             PHSP;
  Enddecay

  Decay eta
  1.000  gamma gamma             PHSP;
  Enddecay

  End
DECAYCARD

exMC_eta_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Ds_eta_enu"
  config.related_dataset = ds_4180
  config.events          = 500_000
  config.decay_card      = decay_card_eta_enu
  config.cross_section   = :default
end

exMC_etap_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Ds_etap_enu"
  config.related_dataset = ds_4180
  config.events          = 500_000
  config.decay_card      = decay_card_etap_enu
  config.cross_section   = :default
end

# ============================================================================
# Algorithm 1: Ds+ -> eta e+ nu_e
# ============================================================================

alg_eta_enu = TagAnalysis.new("DsTagEtaEnu")
alg_eta_enu.set_header(["DsTagEtaEnuAlg/DsTagEtaEnu.h"])
           .set_constant("ECMS" => [:double, 4.178])

alg_eta_enu.tag_side(:Ds) { |t|
  t.modes :all
  t.charm(-1)
}

alg_eta_enu.signal_side { |s|
  s.photons 2
  s.charged(ep: 1, at_least: true)
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
}

alg_eta_enu.fit { |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
}

alg_eta_enu
  .note(:ds_production, "Ds+ mesons produced mainly via e+e- -> Ds+ Ds*- + c.c.; Ds*- -> gamma(pi0) Ds-")
  .note(:st_method, "ST Ds- reconstructed from 14 hadronic decay modes; M_BC in (2.010,2.073) GeV/c^2; select candidate with M_rec closest to Ds*+ mass")
  .note(:st_yield, "Total ST yield: 395142 +/- 1923 from fits to M_tag spectra")
  .note(:gamma_pi0_transition, "Photon or pi0 from Ds*+ transition selected by minimizing |Delta_E|; Delta_E window (-0.04, 0.04) GeV")
  .note(:mm2_fit, "Signal yield from unbinned ML fit to MM^2 distribution; MM^2 = missing mass squared of neutrino")
  .note(:etagamma_extra, "E_extra_max < 0.3 GeV and no extra charged tracks to suppress hadronic backgrounds")
  .note(:fsr_recovery, "Bremsstrahlung energy partially recovered by adding EMC showers within 10 degrees of e+ direction")
  .note(:eta_subdecays, "eta -> gamma gamma (M(γγ) in [0.50,0.57] GeV/c^2, mass-constrained via kinematic fit); eta -> pi+pi-pi0 also used")
  .note(:form_factor, "f_+^{eta/etap}(0)|V_cs| determined from q^2-dependent differential decay width; modified pole model and series expansion fits")
  .note(:eta_eta_p_mixing, "eta-eta' mixing angle phi_P = (40.1 +/- 2.1 +/- 0.7) degrees determined from BF ratios")
  .note(:conexc_mc, "Open charm processes generated using ConExc; ISR and FSR effects considered")
  .with_decay_card(decay_card_eta_enu)
  .apply

# ============================================================================
# Algorithm 2: Ds+ -> eta' e+ nu_e
# ============================================================================

alg_etap_enu = TagAnalysis.new("DsTagEtapEnu")
alg_etap_enu.set_header(["DsTagEtapEnuAlg/DsTagEtapEnu.h"])
            .set_constant("ECMS" => [:double, 4.178])

alg_etap_enu.tag_side(:Ds) { |t|
  t.modes :all
  t.charm(-1)
}

alg_etap_enu.signal_side { |s|
  s.photons 2
  s.charged(pip: 1, pim: 1, ep: 1, at_least: true)
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
}

alg_etap_enu.fit { |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:pip, :pim, :eta).constrain_to_nominal_mass_of(:etap)
  f.chi2_cut 200
}

alg_etap_enu
  .note(:ds_production, "Ds+ mesons produced mainly via e+e- -> Ds+ Ds*- + c.c.")
  .note(:st_method, "ST Ds- reconstructed from 14 hadronic decay modes; same ST selection as eta e nu channel")
  .note(:etap_subdecays, "eta' -> eta pi+ pi- (eta -> gamma gamma) and eta' -> gamma rho0 (rho0 -> pi+ pi-) sub-decays used")
  .note(:mm2_fit, "Signal yield from simultaneous unbinned ML fit to MM^2 spectra of both etap sub-decay modes")
  .note(:helicity_cut, "cos(theta_hel) in (-0.85, 0.85) for etap->gamma rho0 to suppress Ds+ -> etap pi+ and Ds+ -> phi e+ nu backgrounds")
  .note(:metap_cut, "M(etap e+) < 1.9 GeV/c^2 cut to suppress backgrounds")
  .note(:form_factor, "f_+^{etap}(0)|V_cs| measured: 0.477 +/- 0.049(stat) +/- 0.011(syst) (2-parameter series expansion)")
  .note(:peaking_bkg, "Ds+ -> phi e+ nu_e peaking background modeled separately in MM^2 fit for etap -> gamma rho0 channel")
  .with_decay_card(decay_card_etap_enu)
  .apply

# ============================================================================
# Execute
# ============================================================================

root_eta  = alg_eta_enu.execute_on([ds_4180, incMC_4180, exMC_eta_enu])
root_etap = alg_etap_enu.execute_on([ds_4180, incMC_4180, exMC_etap_enu])