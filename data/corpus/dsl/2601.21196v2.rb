# BESIII analysis (companion paper to 2601.21185v2):
# Precise measurements of D^{0(+)} -> Kbar l+ nu_l (l = e, mu) at sqrt(s)=3.773 GeV.
# 20.3 fb^-1 psi(3770) data, double-tag (DT) technique.
# 4 signal channels:
#   D0  -> K-      e+  nu_e
#   D0  -> K-      mu+ nu_mu
#   D+  -> Kbar0   e+  nu_e   (Kbar0 -> K_S0 -> pi+ pi-)
#   D+  -> Kbar0   mu+ nu_mu
# Six ST tag modes for each of Dbar0 and D- (from Sec. IV Table 1).

### Dataset ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

### Decay cards ###
decay_card_D0_Kenu = <<~DECAYCARD
  Decay D*0
  1.0000  D0                                  PHSP;
  Enddecay

  Decay D0
  1.0000  K-  e+  nu_e                        ISGW2;
  Enddecay

  End
DECAYCARD

decay_card_D0_Kmunu = <<~DECAYCARD
  Decay D*0
  1.0000  D0                                  PHSP;
  Enddecay

  Decay D0
  1.0000  K-  mu+  nu_mu                      ISGW2;
  Enddecay

  End
DECAYCARD

decay_card_Dp_K0enu = <<~DECAYCARD
  Decay D+
  1.0000  anti-K0  e+  nu_e                   ISGW2;
  Enddecay

  Decay anti-K0
  1.0000  K_S0                                PHSP;
  Enddecay

  Decay K_S0
  1.0000  pi+  pi-                            PHSP;
  Enddecay

  End
DECAYCARD

decay_card_Dp_K0munu = <<~DECAYCARD
  Decay D+
  1.0000  anti-K0  mu+  nu_mu                 ISGW2;
  Enddecay

  Decay anti-K0
  1.0000  K_S0                                PHSP;
  Enddecay

  Decay K_S0
  1.0000  pi+  pi-                            PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_D0_Kenu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "D0_Kmenu";  c.related_dataset = psi3770_data
  c.events = 2_000_000; c.decay_card = decay_card_D0_Kenu; c.cross_section = :default
end
exMC_D0_Kmunu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "D0_Kmmunu"; c.related_dataset = psi3770_data
  c.events = 2_000_000; c.decay_card = decay_card_D0_Kmunu; c.cross_section = :default
end
exMC_Dp_K0enu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "Dp_K0enu";  c.related_dataset = psi3770_data
  c.events = 2_000_000; c.decay_card = decay_card_Dp_K0enu; c.cross_section = :default
end
exMC_Dp_K0munu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "Dp_K0munu"; c.related_dataset = psi3770_data
  c.events = 2_000_000; c.decay_card = decay_card_Dp_K0munu; c.cross_section = :default
end
exMC_D0_Kenu.save_to_config(format: :yaml,   file_path: 'exMC_D0_Kmenu')
exMC_D0_Kmunu.save_to_config(format: :yaml,  file_path: 'exMC_D0_Kmmunu')
exMC_Dp_K0enu.save_to_config(format: :yaml,  file_path: 'exMC_Dp_K0enu')
exMC_Dp_K0munu.save_to_config(format: :yaml, file_path: 'exMC_Dp_K0munu')

########################################################################
# Algorithm I : D0 -> K- e+ nu_e (Dbar0 hadronic single tag)
# 6 tag modes: Kpi, Kpipi0, K3pi, Kspipi, Kpipi0pi0, K3pipi0.
########################################################################
alg_D0Kenu = TagAnalysis.new("D0toKmEpNue")
alg_D0Kenu.set_header(["D0toKmEpNueAlg/D0toKmEpNue.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })

alg_D0Kenu.tag_side(:D0) do |t|
  t.modes :D0toKPi,
          :D0toKPiPi0,
          :D0toKPiPiPi,
          :D0toKsPiPi,
          :D0toKPiPi0Pi0,
          :D0toKPiPiPiPi0
  t.charm -1                                # Dbar0 tag
end

alg_D0Kenu.signal_side do |s|
  s.charged(km: 1, ep: 1)                   # K- e+ on the signal side
  s.require_charge 0                         # -1 + 1 = 0
  s.missing :nu_e                            # massless neutrino
  s.min_photon_energy 0.025
end

alg_D0Kenu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_D0Kenu
  .note(:tag_deltaE_windows,
        "Per-tag-mode DeltaE windows (MeV) around fitted peaks: Kpi (-27,+27), " \
        "Kpipi0 (-62,+49), K3pi (-26,+24), Kspipi (-24,+24), Kpipi0pi0 " \
        "(-68,+53), K3pipi0 (-57,+51). Applied on tag side.")
  .note(:tag_mBC_window,
        "Tag M_BC required in (1.859, 1.873) GeV/c^2 for Dbar0 tags.")
  .note(:positron_pid,
        "Positron PID: L_e > 0.001 and L_e > 0.8*(L_e + L_pi + L_K).")
  .note(:extra_signal_side_veto,
        "No additional good charged track on the signal side (N_extra_trk = 0); " \
        "maximum energy of extra photons (E_extra_gamma_max) < 0.25 GeV.")
  .note(:mKl_veto,
        "To veto pi+/l+ misID backgrounds: M(K- e+) < 1.83 GeV/c^2 (applied " \
        "in ROOT stage as a stored-variable cut).")
  .note(:Umiss_signal_extraction,
        "Signal yield from binned max-likelihood fit to Umiss=Emiss-|pmiss|c. " \
        "p_D reconstructed as p_D = -phat_Dbar * sqrt(Ebeam^2/c^2 - m_Dbar^2)c^2.")
  .with_decay_card(decay_card_D0_Kenu)
  .apply

alg_D0Kenu.execute_on([psi3770_data, psi3770_incMC, exMC_D0_Kenu])

########################################################################
# Algorithm II : D0 -> K- mu+ nu_mu
########################################################################
alg_D0Kmunu = TagAnalysis.new("D0toKmMupNumu")
alg_D0Kmunu.set_header(["D0toKmMupNumuAlg/D0toKmMupNumu.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })

alg_D0Kmunu.tag_side(:D0) do |t|
  t.modes :D0toKPi,
          :D0toKPiPi0,
          :D0toKPiPiPi,
          :D0toKsPiPi,
          :D0toKPiPi0Pi0,
          :D0toKPiPiPiPi0
  t.charm -1
end

alg_D0Kmunu.signal_side do |s|
  s.charged(km: 1, mup: 1)                  # K- mu+
  s.require_charge 0
  s.missing :nu_mu
  s.min_photon_energy 0.025
end

alg_D0Kmunu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_D0Kmunu
  .note(:tag_deltaE_windows,
        "Per-tag-mode DeltaE windows (MeV): as Algorithm I.")
  .note(:tag_mBC_window,
        "Tag M_BC in (1.859, 1.873) GeV/c^2 for Dbar0 tags.")
  .note(:muon_pid,
        "Muon PID: L_mu > L_e and L_mu > 0.001; also E_EMC(mu) in " \
        "(0.1, 0.3) GeV to suppress e/pi contamination.")
  .note(:extra_signal_side_veto,
        "N_extra_trk = 0 and E_extra_gamma_max < 0.25 GeV.")
  .note(:mKl_veto,
        "To suppress D -> K pi pi0 with pi->mu misID: M(K- mu+) < 1.56 GeV/c^2.")
  .note(:Umiss_signal_extraction,
        "Umiss fit with peaking background from D -> K pi pi0 modeled " \
        "using dedicated MC-simulated shape.")
  .with_decay_card(decay_card_D0_Kmunu)
  .apply

alg_D0Kmunu.execute_on([psi3770_data, psi3770_incMC, exMC_D0_Kmunu])

########################################################################
# Algorithm III : D+ -> Kbar0 e+ nu_e (Kbar0 -> K_S0 -> pi+ pi-)
# 6 tag modes: K+pi-pi-, Ks pi-, K+pi-pi-pi0, Ks pi- pi0, Ks 3pi, KK pi.
########################################################################
alg_DpK0enu = TagAnalysis.new("DptoK0EpNue")
alg_DpK0enu.set_header(["DptoK0EpNueAlg/DptoK0EpNue.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })

alg_DpK0enu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1                                # D- hadronic tag
end

alg_DpK0enu.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)          # pi+ pi- (from K_S0) + e+
  s.require_charge 1
  s.missing :nu_e
  s.min_photon_energy 0.025
end

alg_DpK0enu.fit do |f|
  f.constrain_four_momentum
  # K_S0 mass constraint on the (pi+, pi-) system from Kbar0 -> Ks
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_DpK0enu
  .note(:tag_deltaE_windows,
        "Per-tag-mode DeltaE windows (MeV): K+pi-pi- (-25,+24), Ks pi- (-25,+26), " \
        "K+pi-pi-pi0 (-57,+46), Ks pi- pi0 (-62,+49), Ks 3pi (-28,+27), " \
        "KK pi (-24,+23).")
  .note(:tag_mBC_window,
        "Tag M_BC in (1.863, 1.877) GeV/c^2 for D- tags.")
  .note(:Ks_reconstruction,
        "Kbar0 via K_S0 -> pi+ pi-: 2 oppositely charged tracks (no PID), " \
        "|cos theta|<0.93, |Vz|<20 cm, vertex chi^2<100, decay length > 2*sigma, " \
        "M(pi+pi-) in (0.487, 0.511) GeV/c^2.")
  .note(:positron_pid,
        "Positron PID: L_e > 0.001 and L_e > 0.8*(L_e + L_pi + L_K).")
  .note(:extra_signal_side_veto,
        "N_extra_trk = 0 and E_extra_gamma_max < 0.25 GeV.")
  .note(:mKl_veto,
        "M(Kbar0 e+) < 1.84 GeV/c^2 to suppress pi+/e+ misID backgrounds.")
  .with_decay_card(decay_card_Dp_K0enu)
  .apply

alg_DpK0enu.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_K0enu])

########################################################################
# Algorithm IV : D+ -> Kbar0 mu+ nu_mu
########################################################################
alg_DpK0munu = TagAnalysis.new("DptoK0MupNumu")
alg_DpK0munu.set_header(["DptoK0MupNumuAlg/DptoK0MupNumu.h"])
            .set_constant({ "ECMS" => [:double, 3.773] })

alg_DpK0munu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

alg_DpK0munu.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)
  s.require_charge 1
  s.missing :nu_mu
  s.min_photon_energy 0.025
end

alg_DpK0munu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_DpK0munu
  .note(:tag_deltaE_windows,
        "Per-tag-mode DeltaE windows (MeV): as Algorithm III.")
  .note(:tag_mBC_window,
        "Tag M_BC in (1.863, 1.877) GeV/c^2 for D- tags.")
  .note(:Ks_reconstruction,
        "Kbar0 -> K_S0 -> pi+ pi- with standard BESIII K_S0 selection.")
  .note(:muon_pid,
        "Muon PID: L_mu > L_e, L_mu > 0.001, E_EMC(mu) in (0.1, 0.3) GeV.")
  .note(:extra_signal_side_veto,
        "N_extra_trk = 0 and E_extra_gamma_max < 0.25 GeV.")
  .note(:mKl_veto,
        "M(Kbar0 mu+) < 1.59 GeV/c^2 to suppress D -> Ks pi pi0 backgrounds.")
  .note(:Umiss_signal_extraction,
        "Umiss fit with peaking background from D -> Kbar pi pi0.")
  .with_decay_card(decay_card_Dp_K0munu)
  .apply

alg_DpK0munu.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_K0munu])
