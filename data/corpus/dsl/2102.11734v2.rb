# =============================================================================
# BESIII arXiv:2102.11734v2
# Measurement of absolute branching fractions for purely leptonic Ds+ decays
# at sqrt(s) = 4.178, 4.190, 4.200, 4.210, 4.220, 4.226 GeV (6.32 fb^-1 total).
#
# Method: Double-tag (DT) technique.  The Ds- is reconstructed in 13 hadronic
# decay modes (ST, single-tag sample) and the purely leptonic decay
# Ds+ -> mu+ nu_mu (or Ds+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau) is
# searched for in the system recoiling against the tag.  The neutrino is
# undetected; its kinematics is inferred from four-momentum conservation
# through the missing-mass-squared (MM^2) distribution, which peaks at zero
# for the muonic signal and at a higher value for the tauonic signal
# (two neutrinos are missing).
#
# The absolute branching fraction follows from
#   B(Ds+ -> l+ nu_l) = N_DT / (N_ST_tot x epsilon_DT),
# where epsilon_DT is the DT efficiency corrected for tag-side bias and the
# ST yields are obtained from fits to the M_BC distributions of the 13
# hadronic tag modes.
#
# Charge-conjugate processes are included throughout.
#
# Two independent signal channels are implemented as two separate TagAnalysis
# blocks (Rule T1: different charged multisets).
# =============================================================================

### Dataset description (6 centre-of-mass energies, 6.32 fb^-1 total) ###
data_4180  = DatasetManager.real_data.find("703_4180")    # sqrt(s) = 4.178 GeV
data_4190  = DatasetManager.real_data.find("703_4190")    # sqrt(s) = 4.190 GeV
data_4200  = DatasetManager.real_data.find("703_4200")    # sqrt(s) = 4.200 GeV
data_4210  = DatasetManager.real_data.find("703_4210")    # sqrt(s) = 4.210 GeV
data_4220  = DatasetManager.real_data.find("703_4220")    # sqrt(s) = 4.220 GeV
data_4230  = DatasetManager.real_data.find("703_4230")    # sqrt(s) = 4.226 GeV

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180") # generic inclusive MC (KKMC + EvtGen)
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

data_points = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]
inc_mcs     = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]
datasets    = data_points + inc_mcs

### Decay cards for exclusive signal MC ###
# Signal MC for the muonic channel: Ds+ -> mu+ nu_mu.
# The Ds- decays inclusively over the tag modes (anything).
# The production is through e+e- -> Ds+ Ds- at the psi(4260); the Ds*
# radiative photon is not modelled here (see dsstar_radiative_photon note).
decay_card_mu = <<~DECAYCARD
  Decay psi(4260)
  1.0000  D_s+  D_s-  PHSP;
  Enddecay

  Decay D_s+
  1.0000  mu+  nu_mu  PHSP;
  Enddecay

  Decay D_s-
  1.0000  anything  PHSP;
  Enddecay

  End
DECAYCARD

# Signal MC for the tauonic channel: Ds+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau.
# The signal pi+ is the only charged track on the signal side; both neutrinos
# (nu_tau and anti-nu_tau) are undetected.
decay_card_tau = <<~DECAYCARD
  Decay psi(4260)
  1.0000  D_s+  D_s-  PHSP;
  Enddecay

  Decay D_s+
  1.0000  tau+  nu_tau  PHSP;
  Enddecay

  Decay tau+
  1.0000  pi+  anti-nu_tau  PHSP;
  Enddecay

  Decay D_s-
  1.0000  anything  PHSP;
  Enddecay

  End
DECAYCARD

exMC_mu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "DsLeptonicMuSig"
  config.events        = 500_000
  config.decay_card    = decay_card_mu
  config.cross_section = :default
end
exMC_mu.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_tau = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "DsLeptonicTauSig"
  config.events        = 500_000
  config.decay_card    = decay_card_tau
  config.cross_section = :default
end
exMC_tau.save_to_config(format: :yaml, file_path: 'temp_for_test')

########################################################################
# Algorithm I : Ds+ -> mu+ nu_mu  (muonic channel)
########################################################################
alg_mu = TagAnalysis.new("DsLeptonicMu")
alg_mu.set_header(["DsLeptonicMuAlg/DsLeptonicMu.h"])
      .set_constant({ "ECMS" => [:double, 4.178] })
      .note(:tag_side_selection,
            "The Ds- is reconstructed in 13 hadronic decay modes using the " \
            "pre-stored tag collection.  Charged tracks require " \
            "|cos(theta)| < 0.93 with respect to the beam direction and a " \
            "distance of closest approach to the interaction point within " \
            "10 cm along the beam (V_z) and 1 cm in the perpendicular plane " \
            "(V_xy); tracks from K_S0 and Lambda decays are exempt.  " \
            "Pion/kaon separation combines the MDC dE/dx and TOF information " \
            "into combined confidence levels: a track is identified as a kaon " \
            "when CL_K > CL_pi and as a pion when CL_pi > CL_K.  Photon " \
            "candidates are isolated EMC clusters with E > 25 MeV in the " \
            "barrel (|cos(theta)| < 0.80) and E > 50 MeV in the endcap " \
            "(0.86 < |cos(theta)| < 0.92), with the EMC time in (0, 700) ns " \
            "and more than 10 deg from any charged track.  pi0 candidates " \
            "are photon pairs with 0.115 < M(gamma gamma) < 0.150 GeV/c^2, " \
            "with a 1C kinematic fit constraining the invariant mass to the " \
            "nominal pi0 mass (chi2 < 20).  K_S0 candidates are reconstructed " \
            "from oppositely charged track pairs with a secondary vertex fit, " \
            "positive decay length, and 0.485 < M(pi+ pi-) < 0.510 GeV/c^2.")
      .note(:tag_mbc_deltae,
            "Tag candidates are identified via the beam-constrained mass " \
            "M_BC = sqrt(E_beam^2/c^4 - |p_Ds-|^2/c^2) and the energy " \
            "difference Delta E = E_beam - E_Ds-.  The ST yields are obtained " \
            "from unbinned maximum-likelihood fits to the M_BC distributions " \
            "in the signal region, keeping |Delta E| within approximately " \
            "+-3 sigma around the Delta E peak (mode-dependent windows).  " \
            "The total ST yield summed over the 13 modes is N_ST_tot; the " \
            "per-mode yields and the M_BC signal region (M_BC > 2.05 GeV/c^2) " \
            "are stored and applied at the ROOT stage.")
      .note(:dsstar_radiative_photon,
            "At these centre-of-mass energies (4.178-4.226 GeV) the dominant " \
            "production is e+e- -> Ds*+ Ds- followed by Ds*+ -> gamma Ds+.  " \
            "The radiative photon is selected by requiring its energy in the " \
            "Ds*+ rest frame to lie within 119-149 MeV (the expected value " \
            "from the Ds*+ - Ds+ mass difference).  This photon selection is " \
            "handled at the ROOT analysis stage and is not expressed in the " \
            "DSL; the candidate with photon energy closest to the expected " \
            "value is retained.")
      .note(:muon_pid,
            "The signal muon identification uses probabilities L'_mu, L'_e, " \
            "and L'_K computed from the MDC dE/dx, the TOF and the EMC " \
            "deposited energy.  A muon candidate must satisfy L'_mu > 0.001, " \
            "L'_mu > L'_e and L'_mu > L'_K.  The muon must carry positive " \
            "charge (+1) opposite to the tagged Ds-.  Muons from punch-through " \
            "and decay-in-flight are suppressed by the combined PID requirement.")
      .note(:signal_side_selection,
            "The signal side must contain exactly one positively charged track " \
            "(the muon) and no extra photon beyond the Ds* radiative photon " \
            "(E_extra_gamma_max < 0.3 GeV).  Extra charged tracks are vetoed " \
            "to reject hadronic Ds+ decays.  When multiple tag candidates " \
            "exist the one with the smallest |Delta E| is kept.")
      .note(:background_veto,
            "The dominant background is Ds+ -> tau+ nu_tau (tau+ -> pi+ " \
            "anti-nu_tau), which produces a broad MM^2 distribution.  This " \
            "peaking background is constrained in the simultaneous fit of " \
            "both channels.  Hadronic Ds+ decays with missing neutral " \
            "particles are suppressed by the single-track requirement.")
      .note(:mm2_signal_region,
            "The missing-mass-squared is defined as " \
            "MM^2 = (E_beam - E_l+)^2/c^4 - |p_Ds+ - p_l+|^2/c^2, where " \
            "p_Ds+ = -p_hat_tag sqrt(E_beam^2/c^2 - m^2_Ds- c^2) using the " \
            "ST Ds- direction and the nominal Ds- mass.  For the muonic " \
            "channel the signal peaks at MM^2 ~ 0 (M_nu_mu = 0).  The yield " \
            "is extracted from an unbinned maximum-likelihood fit to the MM^2 " \
            "distribution: a Crystal Ball function for the signal and an " \
            "ARGUS function for the combinatorial background, together with " \
            "the Ds+ -> tau+ nu_tau peaking background whose shape is fixed " \
            "from signal MC.  Applied at the ROOT stage.")
      .note(:tag_mode_bias_correction,
            "The DT efficiency epsilon_DT is the ratio of events that survive " \
            "both the ST and the signal-side selections relative to the ST " \
            "efficiency, weighted by the per-mode ST yields in data.  " \
            "Efficiencies are evaluated with signal MC from the six energy " \
            "points and a weighted average is used.")
      .note(:absolute_branching_fraction,
            "B(Ds+ -> mu+ nu_mu) = N_DT_mu / (N_ST_tot x epsilon_DT_mu).  " \
            "The result, averaged over the six centre-of-mass energies, is " \
            "reported in Table I of the paper; the energy-dependent cross " \
            "section is also extracted.")
      .note(:multi_energy,
            "The analysis combines six data samples at sqrt(s) = 4.178, 4.190, " \
            "4.200, 4.210, 4.220, 4.226 GeV (6.32 fb^-1 total).  The ECMS " \
            "constant is set to 4.178 GeV as a nominal fallback; the actual " \
            "per-run beam energy and boost are read from MeasuredEcmsSvc at " \
            "runtime.  The luminosity-weighted averaging and the Born " \
            "cross-section fit are performed at the ROOT stage.")

# Tag side: Ds+ reconstructed in 13 hadronic modes.
# charm 1 selects the Ds+; the signal side then targets the recoiling Ds-.
alg_mu.tag_side(:Ds) do |t|
  t.modes :all
  t.charm 1
end

# Signal side: the mu+ plus the undetected neutrino.
alg_mu.signal_side do |s|
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

# Kinematic fit: 4-momentum conservation including the missing neutrino.
alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu.with_decay_card(decay_card_mu).apply
alg_mu.execute_on(datasets + [exMC_mu])

########################################################################
# Algorithm II : Ds+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau  (tauonic channel)
########################################################################
alg_tau = TagAnalysis.new("DsLeptonicTau")
alg_tau.set_header(["DsLeptonicTauAlg/DsLeptonicTau.h"])
       .set_constant({ "ECMS" => [:double, 4.178] })
       .note(:tag_side_selection,
             "Same Ds- tag selection as the muonic channel: 13 hadronic " \
             "modes with the same charged-track, photon, pi0 and K_S0 " \
             "criteria.  The ST yield N_ST_tot is common to both channels.")
       .note(:dsstar_radiative_photon,
             "The Ds*+ -> gamma Ds+ radiative photon is selected within " \
             "119-149 MeV in the Ds*+ rest frame, same as the muonic channel.  " \
             "Handled at the ROOT stage.")
       .note(:tau_decay_chain,
             "The decay chain is Ds+ -> tau+ nu_tau followed by " \
             "tau+ -> pi+ anti-nu_tau.  The DSL supports a single missing " \
             "particle on the signal side; :nu_tau is declared as the " \
             "primary neutrino from the Ds+ decay, while the anti-nu_tau " \
             "from the tau+ decay is not explicitly listed.  Its " \
             "four-momentum is absorbed into the missing momentum in the " \
             "four-momentum-conservation constraint, so the reconstructed " \
             "MM^2 observable includes both neutrinos.  Consequently the " \
             "MM^2 signal peak is broader and centred at a positive value " \
             "near m_tau^2 - m_pi^2 rather than at zero.")
       .note(:pion_pid,
             "The signal pion from tau+ -> pi+ anti-nu_tau must satisfy " \
             "pion identification: combined dE/dx and TOF likelihoods with " \
             "CL_pi > Cl_K and CL_pi > CL_p.  The pion must carry positive " \
             "charge (+1) opposite to the tagged Ds-.  Muon veto is applied " \
             "to suppress the Ds+ -> mu+ nu_mu feed-across: the pion must " \
             "fail the muon PID criteria (L'_mu < L'_pi).")
       .note(:signal_side_selection,
             "The signal side must contain exactly one positively charged " \
             "track (the pion from the tau decay) and no extra photon beyond " \
             "the Ds* radiative photon (E_extra_gamma_max < 0.3 GeV).  Extra " \
             "charged tracks are vetoed to reject hadronic Ds+ decays.")
       .note(:background_veto,
             "The dominant background is Ds+ -> mu+ nu_mu, which is " \
             "suppressed by the muon PID veto on the signal track.  Hadronic " \
             "Ds+ decays are rejected by the single-track requirement.  The " \
             "combinatorial background is modelled with an ARGUS function.")
       .note(:mm2_signal_region,
             "The MM^2 distribution for the tauonic channel peaks at " \
             "approximately +1.4 GeV^2/c^4.  The signal yield is extracted " \
             "from a simultaneous unbinned maximum-likelihood fit of the MM^2 " \
             "distributions of both the muonic and tauonic channels, together " \
             "with the Ds+ -> mu+ nu_mu cross-feed whose shape and " \
             "normalisation are constrained by the muonic-channel fit.  " \
             "Applied at the ROOT stage.")
       .note(:absolute_branching_fraction,
             "B(Ds+ -> tau+ nu_tau) = N_DT_tau / (N_ST_tot x epsilon_DT_tau), " \
             "where epsilon_DT_tau includes the tau+ -> pi+ anti-nu_tau " \
             "branching fraction taken from the PDG.  The result, averaged " \
             "over the six centre-of-mass energies, is reported together with " \
             "the muonic-channel result in Table I of the paper.")
       .note(:multi_energy,
             "Same six energy points as the muonic channel (4.178-4.226 GeV), " \
             "combined in the simultaneous fit of both leptonic channels.  " \
             "The ECMS constant is a nominal fallback; the per-run beam " \
             "energy is read from the database at runtime.")

# Tag side: same as the muonic channel.
alg_tau.tag_side(:Ds) do |t|
  t.modes :all
  t.charm 1
end

# Signal side: pi+ from tau decay plus the undetected neutrino.
# The anti-nu_tau from the tau+ decay is not declared — see tau_decay_chain note.
alg_tau.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
  s.missing :nu_tau
end

# Kinematic fit: 4-momentum conservation including the missing nu_tau.
# The anti-nu_tau momentum is absorbed into the missing four-momentum.
alg_tau.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_tau.with_decay_card(decay_card_tau).apply
alg_tau.execute_on(datasets + [exMC_tau])