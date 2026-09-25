# BESIII: ψ(3686) → pbar K+ Σ0 + c.c., partial-wave analysis and
# observation of a new excited Σ state Σ(2330).
# Reconstruction chain: Σ0 → Λ γ, Λ → p π-
# Final state: pbar K+ p π- γ (with charge conjugate)

### Dataset ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

### Decay card ###
decay_card = <<~DECAYCARD
  Decay psi(2S)
  0.5000  anti-p-  K+  Sigma0                  PHSP;
  0.5000  p+       K-  anti-Sigma0             PHSP;
  Enddecay

  Decay Sigma0
  1.0000  gamma  Lambda0                       PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000  gamma  anti-Lambda0                  PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                              HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  anti-p-  pi+                         HypWK;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_pbar_Kp_Sigma0"
  c.related_dataset = psip_data
  c.events          = 500000
  c.decay_card      = decay_card
  c.cross_section   = :default
end

### Algorithm ###
alg = Algorithm.new("PsipPbarKpSigma0")
alg.set_header(["PsipPbarKpSigma0Alg/PsipPbarKpSigma0.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"     # p and K+ (or pbar+K- in c.c.)
      nChrn     ">=2"
      nTot      "==4"
      nNet      "==0"
    }
    .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track    10.0
      nGam              ">=1"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion, :proton]
      identify :pion,   against: [:kaon, :proton]
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    # Nominal 4C kinematic fit under e+e- -> pbar K+ Λ γ hypothesis
    .kinematic_fit([:prm, :kp, :Lambda, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 45
    }
    # Competing hypothesis: no radiative photon (ψ(3686) → pbar K+ Λ)
    .kinematic_fit([:prm, :kp, :Lambda]) {
      constrain_four_momentum
    }
    # Competing hypothesis: two extra photons (ψ(3686) → pbar K+ Λ γ γ)
    .kinematic_fit([:prm, :kp, :Lambda, :gamma, :gamma]) {
      constrain_four_momentum
    }

alg.note(:lambda_mass_window,
         "Λ candidates required to have |M(p π-) - m_Λ| < 7.5 MeV/c^2 " \
         "after common-vertex fit of p π-.")
   .note(:chi2_ordering_veto,
         "χ²_4C for the signal hypothesis (pbar K+ Λ γ) required to be " \
         "smaller than both χ²_4C(pbar K+ Λ) and χ²_4C(pbar K+ Λ γ γ) to " \
         "suppress ψ(3686) → pbar K+ Λ and pbar K+ Λ γ γ backgrounds.")
   .note(:chi_cJ_background_veto,
         "Recoil mass of the photon required to exceed the nominal χ_c2 " \
         "mass by at least 15 MeV/c^2 to suppress ψ(3686) → γ χ_cJ, " \
         "χ_cJ → pbar K+ Λ backgrounds.")
   .note(:helix_correction,
         "Helix-parameter correction applied to charged tracks before the " \
         "4C kinematic fit.")

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([psip_data, psip_incMC, exMC_signal])
