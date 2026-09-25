# BESIII amplitude analysis and BF measurement of J/psi -> gamma eta pi0
# Data: 10087M J/psi. Signal reconstructed from 5 photons, no charged tracks.

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DC
  Decay J/psi
  1.0000 gamma eta pi0 PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DC

exmc_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Jpsi_gammaEtaPi0"
  c.related_dataset = jpsi_data
  c.events          = 1_000_000
  c.decay_card      = decay_card
  c.cross_section   = :default
end

alg = Algorithm.new("JpsiGammaEtaPi0")
alg.set_header(["JpsiGammaEtaPi0Alg/JpsiGammaEtaPi0.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz 10.0
      Vr 1.0
      nChrp "==0"
      nChrn "==0"
      nNet "==0"
    }
   .select_photon {
      tdc_emc_start   0
      tdc_emc_end     14
      angle_to_track  10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=5"
    }
    # 6C kinematic fit: 4-momentum conservation + pi0 and eta mass constraints
    # (best 5-photon combination with smallest chi2)
   .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 200
    }

alg.note(:min_delta_pi0_sq,
         "require min(Delta_pi0^2) = min[(m(g,g)-m_pi0)^2 + (m(g',g'')-m_pi0)^2] > 0.05 GeV^2/c^4 " \
         "over all 4-photon combinations to suppress wrong pi0/eta pairings")
   .note(:omega_veto_pi0geta,
         "veto |m(pi0 gamma_eta) - m_omega| < 65 MeV/c^2 to suppress omega -> gamma pi0")
   .note(:omega_veto_gammageta,
         "veto |m(gamma gamma_eta) - m_omega| < 65 MeV/c^2")
   .note(:eta_veto_gammagammapi0,
         "veto |m(gamma gamma_pi0) - m_eta| < 39 MeV/c^2 to suppress J/psi -> gamma eta eta")
   .note(:pi0_eta_selection,
         "select the 5-photon assignment whose two-photon masses best match nominal pi0 (5 MeV/c^2) " \
         "and eta (9 MeV/c^2) resolutions")
   .note(:tight_4c_chi2,
         "the paper uses the 5-photon combination with the smallest chi2_4C, requiring chi2_4C < 25; " \
         "here the loose default chi2_cut=200 is applied in BOSS, the tight cut is applied in ROOT")

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([jpsi_data, jpsi_incMC, exmc_signal])
