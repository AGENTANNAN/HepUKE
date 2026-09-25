# =============================================================================
# BESIII: psi' -> pbar K+ Sigma0 (and c.c.) & chi_cJ -> pbar K+ Lambda (+c.c.)
# arXiv:1211.5631v2
# =============================================================================

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# --- Decay cards ---
decay_card_pKSigma0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 anti-p- K+ Sigma0 PHSP;
  Enddecay

  Decay Sigma0
  1.0000 gamma Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_pKLambda_chic0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.0000 anti-p- K+ Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_pKLambda_chic1 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.0000 anti-p- K+ Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_pKLambda_chic2 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c2
  1.0000 anti-p- K+ Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Exclusive MC ---
exMC_pKSigma0     = DatasetManager.create_exclusive_mc { |c| c.sample_name="psip_pbKSigma0";     c.related_dataset=psip_data; c.events=200000; c.decay_card=decay_card_pKSigma0;     c.cross_section=:default }
exMC_pKLambda_c0  = DatasetManager.create_exclusive_mc { |c| c.sample_name="psip_gamma_chic0_pbKL"; c.related_dataset=psip_data; c.events=200000; c.decay_card=decay_card_pKLambda_chic0; c.cross_section=:default }
exMC_pKLambda_c1  = DatasetManager.create_exclusive_mc { |c| c.sample_name="psip_gamma_chic1_pbKL"; c.related_dataset=psip_data; c.events=200000; c.decay_card=decay_card_pKLambda_chic1; c.cross_section=:default }
exMC_pKLambda_c2  = DatasetManager.create_exclusive_mc { |c| c.sample_name="psip_gamma_chic2_pbKL"; c.related_dataset=psip_data; c.events=200000; c.decay_card=decay_card_pKLambda_chic2; c.cross_section=:default }

# =============================================================================
# Common event selection: psi' -> gamma pbar K+ Lambda (Lambda -> p pi-)
# For both psi' -> pbar K+ Sigma0 (Sigma0 -> gamma Lambda) and psi' -> gamma chi_cJ.
# =============================================================================
def build_sel
  sel = Selection.new
  sel.select_track {
       cos_theta 0.93
       Vz  30.0   # |Vz| < 30 cm
       Vr  15.0   # |Vr| < 15 cm
       nChrp ">=2"
       nChrn ">=2"
     }
     .select_photon {
       nGam ">=1"
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       tdc_emc_start 0
       tdc_emc_end 14
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       identify :kaon,   against: [:pion, :proton]
       identify :pion,   against: [:kaon, :proton]
       nprp "==1"
       nprm "==1"
       nkp  "==1"
       npim "==1"
     }
     .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:gamma, :prm, :kp, :Lambda]) {
       nominal
       vertex_fit([1, 2])   # constrain pbar (index 1) and K+ (index 2) to a common vertex
       constrain_four_momentum
       chi2_cut 200
     }
  sel
end

# --- ALG A: psi' -> pbar K+ Sigma0 (Sigma0 -> gamma Lambda) ---
alg_S0 = Algorithm.new("PsipPbKpSigma0")
alg_S0.set_header(["PsipPbKpSigma0Alg/PsipPbKpSigma0.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
alg_S0.note(:pk_common_vertex_tight,
            "pbar+K+ vertex-fit tracks require |Vz|<10 cm and Vr<1 cm (tighter than the "\
            "MDC-only cut) for the gamma pbar K+ Lambda topology")
      .note(:lambda_mass_window,
            "|M(p pi-) - M_Lambda_PDG| < 7 MeV/c^2 (post-4C-fit invariant mass)")
      .note(:sigma0_mass_window,
            "|M(gamma Lambda) - M_Sigma0_PDG| < 15 MeV/c^2 to select Sigma0")
alg_S0.with_decay_card(decay_card_pKSigma0).apply(build_sel)
alg_S0.execute_on([psip_data, psip_incMC, exMC_pKSigma0])

# --- ALG B: psi' -> gamma chi_cJ -> gamma pbar K+ Lambda (J=0,1,2 share the selection)
alg_cJ = Algorithm.new("PsipGammaChicJpbKL")
alg_cJ.set_header(["PsipGammaChicJpbKLAlg/PsipGammaChicJpbKL.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
alg_cJ.note(:sigma0_veto,
            "veto |M(gamma Lambda) - M_Sigma0_PDG| < 15 MeV/c^2 to remove "\
            "psi' -> pbar K+ Sigma0 events")
      .note(:lambda_mass_window,
            "|M(p pi-) - M_Lambda_PDG| < 7 MeV/c^2")
      .note(:pk_common_vertex_tight,
            "pbar+K+ vertex-fit tracks require |Vz|<10 cm and Vr<1 cm")
alg_cJ.with_decay_card(decay_card_pKLambda_chic0).apply(build_sel)
alg_cJ.execute_on([psip_data, psip_incMC, exMC_pKLambda_c0, exMC_pKLambda_c1, exMC_pKLambda_c2])
