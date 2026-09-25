# BESIII: First observation of η_c → Ξ0 Ξbar0
# Decay chain: J/ψ → γ η_c, η_c → Ξ0 Ξbar0, Ξ0 → Λ π0, Ξbar0 → Λbar π0,
#              Λ → p π-, Λbar → pbar π+, π0 → γγ
# Final state: γ p pbar π+ π- γγγγ (5 photons + 4 charged tracks)

### Dataset ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

### Decay card ###
decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000  gamma  eta_c                        JPE;
  Enddecay

  Decay eta_c
  1.0000  Xi0  anti-Xi0                       PHSP;
  Enddecay

  Decay Xi0
  1.0000  Lambda0  pi0                        PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000  anti-Lambda0  pi0                   PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                             HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  anti-p-  pi+                        HypWK;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma                        PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Jpsi_gamma_etac_Xi0Xi0bar"
  c.related_dataset = jpsi_data
  c.events          = 500000
  c.decay_card      = decay_card
  c.cross_section   = :default
end

### Algorithm ###
alg = Algorithm.new("EtacXi0Xi0bar")
alg.set_header(["EtacXi0Xi0barAlg/EtacXi0Xi0bar.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      nChrp   ">=2"       # p and pi+
      nChrn   ">=2"       # pbar and pi-
      nTot    ">=4"
      nNet    "==0"
    }
    .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track    10.0
      nGam              ">=5"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :pion,   against: [:kaon]
    }
    .select_isolated_photon {
      angle_to_prm_track 20.0
      nGam              ">=5"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0 ">=2"
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    # Nominal 4C kinematic fit on γ Λ Λbar π0 π0 (loose χ² cut)
    .kinematic_fit([:gamma, :Lambda, :Lambda_bar, :pi0, :pi0]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }
    # Competing hypothesis: Λ Λbar π0 π0 (no radiative photon) — used for
    # J/ψ → Ξ0 Ξbar0 background suppression (χ² stored, cut applied in ROOT)
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) {
      constrain_four_momentum
    }

alg.note(:lambda_mass_window,
         "Λ (Λbar) reconstructed from p π- (pbar π+); common-vertex fit " \
         "χ² < 200, mass window M(pπ-) ∈ (1.111, 1.120) GeV/c^2.")
   .note(:proton_pion_split,
         "Charged tracks with momentum > 0.247 GeV/c treated as protons, " \
         "< 0.247 GeV/c as pions (Λ decay-kinematics-based separation).")
   .note(:pi0_gamma_gamma_window,
         "π0 candidates formed from γγ pairs with M(γγ) ∈ (0.098, 0.165) " \
         "GeV/c^2 before the 1C mass-constraint kinematic fit.")
   .note(:xi0_mass_window,
         "Best Λπ0 / Λbarπ0 pairing chosen by minimizing " \
         "(m_{π0_1 Λ} - M_{Ξ0})^2 + (m_{π0_2 Λbar} - M_{Ξbar0})^2; " \
         "Ξ0/Ξbar0 mass window (1.296, 1.331) GeV/c^2 (~ ±3σ).")
   .note(:helix_correction,
         "Helix-parameter correction applied to all charged tracks before " \
         "the 4C kinematic fit; systematic estimated with/without correction.")

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
