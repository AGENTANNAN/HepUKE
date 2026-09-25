# BESIII analysis: Measurements of D_s^+ -> l^+ nu_l (l=mu,tau) via e^+e^- -> D_s*+ D_s*-
# arXiv: 2407.11727v3
# Multi-energy double-tag (DT) analysis using 10.64 fb^-1 at sqrt(s) = 4.237-4.699 GeV

### Datasets ###
data_4237  = DatasetManager.real_data.find("703_4237")
data_4246  = DatasetManager.real_data.find("703_4246")
data_4260  = DatasetManager.real_data.find("703_4260")
data_4270  = DatasetManager.real_data.find("703_4270")
data_4280  = DatasetManager.real_data.find("703_4280")
data_4290  = DatasetManager.real_data.find("705_4290")
data_4310  = DatasetManager.real_data.find("705_4310")
data_4315  = DatasetManager.real_data.find("705_4315")
data_4400  = DatasetManager.real_data.find("705_4400")
data_4420  = DatasetManager.real_data.find("703_4420")
data_4440  = DatasetManager.real_data.find("705_4440")
data_4610  = DatasetManager.real_data.find("706_4610")
data_4620  = DatasetManager.real_data.find("706_4620")
data_4640  = DatasetManager.real_data.find("706_4640")
data_4660  = DatasetManager.real_data.find("706_4660")
data_4680  = DatasetManager.real_data.find("706_4680")
data_4700  = DatasetManager.real_data.find("706_4700")

all_energy_points = [data_4237, data_4246, data_4260, data_4270, data_4280,
                     data_4290, data_4310, data_4315, data_4400, data_4420,
                     data_4440, data_4610, data_4620, data_4640, data_4660,
                     data_4680, data_4700]

incMC_4237 = DatasetManager.inclusive_mc.find("703_4237")
incMC_4246 = DatasetManager.inclusive_mc.find("703_4246")
incMC_4260  = DatasetManager.inclusive_mc.find("703_4260")
incMC_4270  = DatasetManager.inclusive_mc.find("703_4270")
incMC_4280  = DatasetManager.inclusive_mc.find("703_4280")
incMC_4290  = DatasetManager.inclusive_mc.find("705_4290")
incMC_4310  = DatasetManager.inclusive_mc.find("705_4310")
incMC_4315  = DatasetManager.inclusive_mc.find("705_4315")
incMC_4400  = DatasetManager.inclusive_mc.find("705_4400")
incMC_4420  = DatasetManager.inclusive_mc.find("703_4420")
incMC_4440  = DatasetManager.inclusive_mc.find("705_4440")
incMC_4610  = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620  = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640  = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660  = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680  = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700  = DatasetManager.inclusive_mc.find("706_4700")

all_incMC = [incMC_4237, incMC_4246, incMC_4260, incMC_4270, incMC_4280,
             incMC_4290, incMC_4310, incMC_4315, incMC_4400, incMC_4420,
             incMC_4440, incMC_4610, incMC_4620, incMC_4640, incMC_4660,
             incMC_4680, incMC_4700]

### Decay cards ###
# Signal process: e+e- -> D_s*+ D_s*-; D_s*+ -> signal, D_s*- -> tag
# Top mother psi(4260) per BESIII KKMC convention

# Decay card for D_s^+ -> mu^+ nu_mu signal
decay_card_ds_munu = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s*- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 mu+ nu_mu SLN;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for D_s^+ -> tau^+ nu_tau, tau^+ -> e^+ nu_e anti-nu_tau
decay_card_ds_taunu_e = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s*- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 tau+ nu_tau SLN;
    Enddecay

    Decay tau+
    1.000 e+ nu_e anti-nu_tau PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for D_s^+ -> tau^+ nu_tau, tau^+ -> mu^+ nu_mu anti-nu_tau
decay_card_ds_taunu_mu = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s*- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 tau+ nu_tau SLN;
    Enddecay

    Decay tau+
    1.000 mu+ nu_mu anti-nu_tau PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for D_s^+ -> tau^+ nu_tau, tau^+ -> pi^+ anti-nu_tau
decay_card_ds_taunu_pi = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s*- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 tau+ nu_tau SLN;
    Enddecay

    Decay tau+
    1.000 pi+ anti-nu_tau PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for D_s^+ -> tau^+ nu_tau, tau^+ -> rho^+ anti-nu_tau, rho^+ -> pi^+ pi0
decay_card_ds_taunu_rho = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s*- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 tau+ nu_tau SLN;
    Enddecay

    Decay tau+
    1.000 rho+ anti-nu_tau PHSP;
    Enddecay

    Decay rho+
    1.000 pi+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC ###
exMC_ds_munu = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "ds_munu_signal"
  config.events        = 1_000_000
  config.decay_card    = decay_card_ds_munu
  config.cross_section = :default
end

exMC_ds_taunu_e = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "ds_taunu_e_signal"
  config.events        = 1_000_000
  config.decay_card    = decay_card_ds_taunu_e
  config.cross_section = :default
end

exMC_ds_taunu_mu = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "ds_taunu_mu_signal"
  config.events        = 1_000_000
  config.decay_card    = decay_card_ds_taunu_mu
  config.cross_section = :default
end

exMC_ds_taunu_pi = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "ds_taunu_pi_signal"
  config.events        = 1_000_000
  config.decay_card    = decay_card_ds_taunu_pi
  config.cross_section = :default
end

exMC_ds_taunu_rho = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "ds_taunu_rho_signal"
  config.events        = 1_000_000
  config.decay_card    = decay_card_ds_taunu_rho
  config.cross_section = :default
end

### Tag-based analysis 1: D_s^+ -> mu^+ nu_mu (mu_a with MUC info) ###
alg_ds_munu = TagAnalysis.new("DsMuNu")
alg_ds_munu.set_header(["DsMuNuAlg/DsMuNu.h"])
  .note(:muon_id_muc, "Muon candidates identified with MUC information (mu_a); depth requirement per
    |cos_theta_mu| and p_mu dependent cuts as in Table III. Signal yield extracted from unbinned
    maximum likelihood fit to M_miss^2 distribution.")
  .note(:tag_deltae, "DeltaE windows per tag mode per Table I. Keep candidate with minimum |DeltaE|
    per mode per charge per event.")
  .note(:tag_mbc_fit, "Single-tag yields from M_BC fit: MC shape convolved with Gaussian (signal),
    ARGUS function (bkg below 4.450 GeV), 2nd-order Chebyshev (bkg above 4.450 GeV).")
  .note(:ntracks_0, "Require N_extra^charge = 0 (no extra charged tracks beyond the signal track).")
  .note(:peaking_veto, "For D_s^- -> pi^+pi^-pi^- and D_s^- -> K^+pi^+pi^- tag modes, veto
    |M(pi^+pi^-) - m(K_S^0)| < 0.03 GeV/c^2 to suppress K_S^0 -> pi^+pi^- background.")

alg_ds_munu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi, :DstoKsK, :DstoKsKPi0,
          :DstoKPiPi, :DstoKsKsPi, :DstoKsKplusPiPi, :DstoKsKminusPiPi,
          :DstoPiEta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam, :DstoPiPi0Eta
  t.charm -1
end

alg_ds_munu.signal_side do |s|
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_ds_munu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_ds_munu.with_decay_card(decay_card_ds_munu).apply
alg_ds_munu.execute_on(all_energy_points + all_incMC + exMC_ds_munu)

### Tag-based analysis 2: D_s^+ -> tau^+ nu_tau (tau^+ -> e^+ nu_e anti-nu_tau) ###
alg_ds_taunu_e = TagAnalysis.new("DsTauENu")
alg_ds_taunu_e.set_header(["DsTauENuAlg/DsTauENu.h"])
  .note(:tau_fit_simultaneous, "Tau BF determined from simultaneous fit to 4 tau decay modes
    (e, mu, pi, rho). E_sum^extra_gamma used for e/mu modes; M_miss^2 used for pi/rho modes.
    Lepton universality constraint R = 9.75 applied in SM-constrained fit.")
  .note(:tag_deltae, "DeltaE windows per tag mode per Table I. Keep candidate with minimum |DeltaE|
    per mode per charge per event.")
  .note(:tag_mbc_fit, "Single-tag yields from M_BC fit: MC shape convolved with Gaussian (signal),
    ARGUS function (bkg below 4.450 GeV), 2nd-order Chebyshev (bkg above 4.450 GeV).")
  .note(:ntracks_0, "Require N_extra^charge = 0 (no extra charged tracks beyond the signal track).")
  .note(:peaking_veto, "For D_s^- -> pi^+pi^-pi^- and D_s^- -> K^+pi^+pi^- tag modes, veto
    |M(pi^+pi^-) - m(K_S^0)| < 0.03 GeV/c^2 to suppress K_S^0 -> pi^+pi^- background.")

alg_ds_taunu_e.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi, :DstoKsK, :DstoKsKPi0,
          :DstoKPiPi, :DstoKsKsPi, :DstoKsKplusPiPi, :DstoKsKminusPiPi,
          :DstoPiEta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam, :DstoPiPi0Eta
  t.charm -1
end

alg_ds_taunu_e.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_tau
end

alg_ds_taunu_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_ds_taunu_e.with_decay_card(decay_card_ds_taunu_e).apply
alg_ds_taunu_e.execute_on(all_energy_points + all_incMC + exMC_ds_taunu_e)

### Tag-based analysis 3: D_s^+ -> tau^+ nu_tau (tau^+ -> mu^+ nu_mu anti-nu_tau) ###
alg_ds_taunu_mu = TagAnalysis.new("DsTauMuNu")
alg_ds_taunu_mu.set_header(["DsTauMuNuAlg/DsTauMuNu.h"])
  .note(:tau_mu_pid, "Muon from tau decay: E_EMC in (0.1, 0.3) GeV; MUC depth requirement
    per |cos_theta| and p_mu as in Table III. Without MUC info (mu_b) for cross-check.")
  .note(:tau_fit_simultaneous, "Tau BF determined from simultaneous fit to 4 tau decay modes.")
  .note(:tag_deltae, "DeltaE windows per tag mode per Table I. Keep candidate with minimum |DeltaE|
    per mode per charge per event.")
  .note(:tag_mbc_fit, "Single-tag yields from M_BC fit.")
  .note(:ntracks_0, "Require N_extra^charge = 0.")
  .note(:peaking_veto, "Veto |M(pi^+pi^-) - m(K_S^0)| < 0.03 GeV/c^2 for pi^+pi^-pi^- and K^+pi^+pi^- tag modes.")

alg_ds_taunu_mu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi, :DstoKsK, :DstoKsKPi0,
          :DstoKPiPi, :DstoKsKsPi, :DstoKsKplusPiPi, :DstoKsKminusPiPi,
          :DstoPiEta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam, :DstoPiPi0Eta
  t.charm -1
end

alg_ds_taunu_mu.signal_side do |s|
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_tau
end

alg_ds_taunu_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_ds_taunu_mu.with_decay_card(decay_card_ds_taunu_mu).apply
alg_ds_taunu_mu.execute_on(all_energy_points + all_incMC + exMC_ds_taunu_mu)

### Tag-based analysis 4: D_s^+ -> tau^+ nu_tau (tau^+ -> pi^+ anti-nu_tau) ###
alg_ds_taunu_pi = TagAnalysis.new("DsTauPiNu")
alg_ds_taunu_pi.set_header(["DsTauPiNuAlg/DsTauPiNu.h"])
  .note(:tau_pi_pid, "Pion from tau decay identified with SimplePIDSvc. Large D_s^+ -> mu^+ nu_mu
    background due to poor mu/pion separation without MUC info. Signal shape: sum of two
    bifurcated-Gaussian functions.")
  .note(:tau_fit_simultaneous, "Tau BF determined from simultaneous fit to 4 tau decay modes.")
  .note(:tag_deltae, "DeltaE windows per tag mode. Keep candidate with minimum |DeltaE| per mode per charge per event.")
  .note(:tag_mbc_fit, "Single-tag yields from M_BC fit.")
  .note(:ntracks_0, "Require N_extra^charge = 0.")
  .note(:peaking_veto, "Veto |M(pi^+pi^-) - m(K_S^0)| < 0.03 GeV/c^2 for pi^+pi^-pi^- and K^+pi^+pi^- tag modes.")

alg_ds_taunu_pi.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi, :DstoKsK, :DstoKsKPi0,
          :DstoKPiPi, :DstoKsKsPi, :DstoKsKplusPiPi, :DstoKsKminusPiPi,
          :DstoPiEta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam, :DstoPiPi0Eta
  t.charm -1
end

alg_ds_taunu_pi.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
  s.missing :nu_tau
end

alg_ds_taunu_pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_ds_taunu_pi.with_decay_card(decay_card_ds_taunu_pi).apply
alg_ds_taunu_pi.execute_on(all_energy_points + all_incMC + exMC_ds_taunu_pi)

### Tag-based analysis 5: D_s^+ -> tau^+ nu_tau (tau^+ -> rho^+ anti-nu_tau, rho^+ -> pi^+ pi0) ###
alg_ds_taunu_rho = TagAnalysis.new("DsTauRhoNu")
alg_ds_taunu_rho.set_header(["DsTauRhoNuAlg/DsTauRhoNu.h"])
  .note(:tau_rho_reco, "rho^+ -> pi^+ pi0, pi0 -> gamma gamma on signal side.
    Signal shape: sum of two bifurcated-Gaussian functions. Main peaking backgrounds:
    D_s^+ -> eta pi^+ pi0 and D_s^+ -> K^0 pi^+ pi0.")
  .note(:tau_fit_simultaneous, "Tau BF determined from simultaneous fit to 4 tau decay modes.")
  .note(:tag_deltae, "DeltaE windows per tag mode. Keep candidate with minimum |DeltaE| per mode per charge per event.")
  .note(:tag_mbc_fit, "Single-tag yields from M_BC fit.")
  .note(:ntracks_0, "Require N_extra^charge = 0.")
  .note(:peaking_veto, "Veto |M(pi^+pi^-) - m(K_S^0)| < 0.03 GeV/c^2 for pi^+pi^-pi^- and K^+pi^+pi^- tag modes.")

alg_ds_taunu_rho.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi, :DstoKsK, :DstoKsKPi0,
          :DstoKPiPi, :DstoKsKsPi, :DstoKsKplusPiPi, :DstoKsKminusPiPi,
          :DstoPiEta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam, :DstoPiPi0Eta
  t.charm -1
end

alg_ds_taunu_rho.signal_side do |s|
  s.photons 2
  s.charged(pip: 1)
  s.require_charge 1
  s.missing :nu_tau
  s.min_photon_angle 10.0
end

alg_ds_taunu_rho.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_ds_taunu_rho.with_decay_card(decay_card_ds_taunu_rho).apply
alg_ds_taunu_rho.execute_on(all_energy_points + all_incMC + exMC_ds_taunu_rho)