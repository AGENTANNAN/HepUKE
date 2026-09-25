### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for psi(3686) -> gamma chi_cJ -> gamma Sigma+ anti-p- K_S0
# Sigma+ -> p+ pi0, pi0 -> gamma gamma, K_S0 -> pi+ pi-
# chi_c0, chi_c1, chi_c2 share the same final state
decay_card_chi_cJ_Sigma = <<~DECAYCARD
    Decay psi(2S)
    0.333 gamma chi_c0 P2GC0;
    0.333 gamma chi_c1 P2GC1;
    0.334 gamma chi_c2 P2GC2;
    Enddecay
    Decay chi_c0
    1.000 Sigma+ anti-p- anti-K_S0 PHSP;
    Enddecay
    Decay chi_c1
    1.000 Sigma+ anti-p- anti-K_S0 PHSP;
    Enddecay
    Decay chi_c2
    1.000 Sigma+ anti-p- anti-K_S0 PHSP;
    Enddecay
    Decay Sigma+
    1.000 p+ pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay anti-K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC sample
exMC_chi_cJ_Sigma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_chi_cJ_to_Sigma_pbar_Ks"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_chi_cJ_Sigma
  config.cross_section = :default
end

### psi(3686) -> gamma chi_cJ -> gamma Sigma+ anti-p- K_S0 ###
# Three chi_cJ states share identical final state (gamma p pbar pi+ pi- gamma gamma)
# Single Algorithm instance serves all three

alg_chi_cJ_Sigma = Algorithm.new("ChiCJToSigmaPbarKs")
alg_chi_cJ_Sigma.set_header(["ChiCJToSigmaPbarKsAlg/ChiCJToSigmaPbarKs.h"])
                 .set_constant({"ECMS" => [:double, 3.686]})
                 .note(:v0_no_vertex_cuts, "No R_xy or V_z requirements on tracks used to form K_S0 and Sigma+ candidates due to their long lifetimes")
                 .note(:background_veto, "Lambda veto: |M(pbar pi+) - m_Lambda| > 6 MeV/c^2; competing hypothesis veto: chi2_4C(gamma gamma gamma) < chi2_4C(gamma gamma) AND chi2_4C(gamma gamma gamma gamma)")
                 .note(:pi0_selection, "pi0 candidate is gamma-gamma pair with invariant mass closest to nominal pi0 mass after 4C fit; third photon is radiative")

sel_chi_cJ_Sigma = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=2"
    nChrn ">=2"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Nominal 4C fit under p pbar pi+ pi- gamma gamma gamma hypothesis
  .kinematic_fit([:prp, :prm, :K_S0, :pi0, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 50
  }
  # Competing hypothesis 4C: p pbar pi+ pi- gamma gamma (no radiative photon) -- stores chi2 for ROOT veto
  .kinematic_fit([:prp, :prm, :K_S0, :pi0, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Competing hypothesis 4C: p pbar pi+ pi- gamma gamma gamma gamma (extra photon) -- stores chi2 for ROOT veto
  .kinematic_fit([:prp, :prm, :K_S0, :pi0, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg_chi_cJ_Sigma.with_decay_card(decay_card_chi_cJ_Sigma).apply(sel_chi_cJ_Sigma)
alg_chi_cJ_Sigma.execute_on([psip_data, psip_incMC, exMC_chi_cJ_Sigma])