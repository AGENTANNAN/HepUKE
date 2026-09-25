### Dataset description ###
# BESIII data at sqrt(s) = 3.773 GeV (psi(3770)), L = 2.93 fb^-1.
# Double-tag (DT) analysis: the D- is reconstructed in six hadronic ST modes and the
# signal D+ -> gamma e+ nu_e is searched for in the remaining tracks and showers.
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process D+ -> gamma e+ nu_e (top mother psi(3770), KKMC)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 gamma e+ nu_e PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the dominant background D+ -> pi0 e+ nu_e, simulated exclusively
decay_card_bkg_pi0enu = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 pi0 e+ nu_e PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DplusToGammaENu"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_pi0enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DplusToPi0ENu"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_bkg_pi0enu
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg = TagAnalysis.new("DplusToGammaENu")
alg.set_header(["DplusToGammaENuAlg/DplusToGammaENu.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: D- (the charge-conjugate partner of the signal D+), reconstructed from the
# pre-stored EvtRecDTag collection in the six hadronic ST modes. Both tag charges are
# scanned (charge conjugate modes always implied); the resolution of K_S0 daughters
# and of pi0 -> gamma gamma is handled by DTagAlg. mBC / deltaE are stored
# unconditionally and windowed in the ROOT analysis.
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,        # K+ pi- pi-
          :DptoKPiPiPi0,     # K+ pi- pi- pi0
          :DptoKsPi,         # K_S0 pi-
          :DptoKsPiPi0,      # K_S0 pi- pi0
          :DptoKsPiPiPi,     # K_S0 pi+ pi- pi-
          :DptoKKPi          # K+ K- pi-
end

# Signal side: exactly one remaining charged track (the electron, charge opposite to the
# tag) and at least one remaining photon (the candidate radiative photon; the highest
# energy photon is used). The undetected neutrino is declared as a massless missing
# particle, which auto-stores m_P4_miss_fit, m_Umiss, m_Umiss2 and m_q2.
alg.signal_side do |s|
  s.charged(ep: 1)          # exactly one leftover track, identified as an electron
  s.photons 1               # >= 1 good remaining shower for the radiative photon
  s.min_photon_angle 10.0
  s.min_photon_energy 0.010 # radiative photon energy > 10 MeV (infrared cutoff)
  s.missing :nu_e           # massless missing neutrino
end

alg.fit do |f|
  f.constrain_four_momentum  # 4C: tag + gamma + e+ + nu_e = ecms_lab
  f.chi2_cut 200
end

alg.note(:signal_mc_model, "Signal MC generated adopting the first-order alpha_s and heavy-quark-mass model of Ref. [11] (Yang & Yang); the minimum energy of the radiative photon is set at 10 MeV to avoid the infrared divergence for soft photons.")
     .note(:tag_side_deltaE, "The ST D- candidates are required to lie in mode-dependent DeltaE regions (MeV): K+pi-pi- [-27,25]; K+pi-pi-pi0 [-62,34]; K_S0 pi- [-25,25]; K_S0 pi- pi0 [-73,41]; K_S0 pi+pi-pi- [-33,30]; K+ K- pi- [-23,20]. These are per-mode windows that the DSL tag_side window cannot express (one window per side), so deltaE is stored by the tag and windowed per mode in the ROOT analysis.")
     .note(:best_candidate_st, "If multiple ST candidates are found in a mode, only the one with the smallest |DeltaE| is retained; the ST M_BC signal region is 1.8628 < M_BC < 1.8788 GeV/c^2.")
     .note(:electron_pid, "The signal-side track is identified as an electron from dE/dx, TOF and EMC with L(e) > 0 and L(e)/(L(e)+L(pi)+L(K)) > 0.8; the lepton PID thresholds in the generated tag code are fixed v1 defaults and are not DSL-tunable.")
     .note(:fsr_recovery, "Energy of neighbouring photons from final-state radiation is added back to the electron candidate: photons with energy > 50 MeV within a 5 degree cone around the electron direction (excluding the radiative photon) are merged into the electron.")
     .note(:pi0_veto, "Background veto: events are rejected if any pair of photons satisfies chi^2 < 20 in the pi0 1C kinematic fit (suppresses D+ -> pi0 e+ nu_e where the two pi0 photons fake the radiative photon).")
     .note(:lateral_moment, "The radiative photon is required to have a lateral moment in (0.0, 0.3) to suppress D+ -> K_L0 e+ nu_e background, where the K_L0 shower shape varies broadly from 0 to 0.85.")
     .note(:bkg_pi0enu_control_sample, "The D+ -> pi0 e+ nu_e background shape is taken from an exclusive MC sample and normalized to an expected yield N_pi0_exp extracted from a data control sample selected with the same ST and electron criteria but requiring a pi0 from two photons with a 1C kinematic fit (chi^2 < 20).")
     .note(:ub_fit, "Signal yield extracted from an extended unbinned maximum-likelihood fit to U_miss = E_miss - |p_miss| c (signal shape from signal MC convoluted with a Gaussian, D+ -> pi0 e+ nu_e shape from exclusive MC, remaining background shape from inclusive MC); performed in the ROOT analysis.")
     .note(:upper_limit, "No significant signal observed; a 90% C.L. upper limit on B(D+ -> gamma e+ nu_e) is obtained with a toy-MC method incorporating statistical and systematic uncertainties (ROOT analysis).")
     .note(:st_common_selection, "Common ST track quality cuts: |cos(theta)| < 0.93, Vr < 1 cm, |Vz| < 10 cm; pi/K separation by PID likelihoods L(pi) > L(K) (pi) and L(K) > L(pi) (K). K_S0: two oppositely charged tracks with |cos(theta)| < 0.93, |Vz| < 20 cm, no Vr and no PID, 0.487 < M(pi+pi-) < 0.511 GeV/c^2, decay length > 2 sigma, momenta after the vertex fit used downstream.")
     .note(:photon_reconstruction, "Photon candidates: EMC showers not associated with charged tracks, E > 25 MeV in the barrel (|cos(theta)| < 0.80) or E > 50 MeV in the endcaps (0.84 < |cos(theta)| < 0.92), shower time within 700 ns of the event start. pi0 -> gamma gamma candidates require 0.115 < M(gamma gamma) < 0.150 GeV/c^2, pairs with both photons in the endcaps are rejected, and a 1C mass-constrained fit is applied.")
     .with_decay_card(decay_card_signal)
     .apply

alg.execute_on([data_3773, incMC_3773, exMC_signal, exMC_pi0enu])
