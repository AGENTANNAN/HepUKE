### Dataset preparation ###
# ψ(3686) at 3.686 GeV
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # ψ(3686) inclusive MC

# ---- Decay cards (EvtGen) ----
# Mode 1: ψ(3686) → γ χ_c2, χ_c2 → K+ K- π0
decay_card_kkpi0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c2
    1.000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2: ψ(3686) → γ χ_c2, χ_c2 → K_S K± π∓ (K_S → π+ π-), both charge configurations
decay_card_kskpi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c2
    0.500 K_S0 K+ pi- PHSP;
    0.500 K_S0 K- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3: ψ(3686) → γ χ_c2, χ_c2 → π+ π- π0
decay_card_pipipi0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c2
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC: 200k events per mode ----
exMC_kkpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic2_kkpi0"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_kkpi0
  config.cross_section   = :default
end

exMC_kskpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic2_kskpi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_kskpi
  config.cross_section   = :default
end

exMC_pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic2_pipipi0"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_pipipi0
  config.cross_section   = :default
end


### Event selection (BOSS) ###

# ---------------- Mode 1: χ_c2 → K+ K- π0 ----------------
alg_name_m1 = "ChiC2GammaKKPi0"
alg_m1 = Algorithm.new(alg_name_m1)
alg_m1.set_header(["#{alg_name_m1}Alg/#{alg_name_m1}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

sel_m1 = Selection.new
  .select_track {                      # charged-track quality cuts + multiplicity
    cos_theta 0.93                     # |cosθ| < 0.93
    Vz 10.0                            # |Vz| < 10 cm
    Vr 1.0                             # Vr < 1 cm
    nChrp "==1"                        # exactly one positive track
    nChrn "==1"                        # exactly one negative track
    nNet  "==0"                        # net charge zero
  }
  .select_photon {                     # photon (EMC shower) selection
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0                # > 10° from any charged track
    energyThreshold_b 0.025            # 25 MeV (barrel)
    energyThreshold_e 0.050            # 50 MeV (endcap)
    nGam ">=3"                         # at least three photons (radiative γ + π0 → γγ)
  }
  .assign({:chrgp => :kp, :chrgn => :km})   # no PID: positive track → K+, negative track → K-
  .kalman_kinematic_fit([:gamma, :gamma]) { # π0 from two photons
    invariant_mass_of(:gamma, :gamma).within(0.125, 0.145)             # |Mγγ - mπ0| < 10 MeV
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                                                         # at least one π0
  }
  .kinematic_fit([:gamma, :kp, :km, :pi0]) {   # 4C fit of γ K+ K- π0
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km).within(0.0, 3.0)  # M(K+K-) in 0-3 GeV
    chi2_cut 80
  }

alg_m1.with_decay_card(decay_card_kkpi0).apply(sel_m1)
root_files_m1 = alg_m1.execute_on([psip_data, psip_incMC, exMC_kkpi0])


# ---------------- Mode 2: χ_c2 → K_S K± π∓ ----------------
alg_name_m2 = "ChiC2GammaKSKPi"
alg_m2 = Algorithm.new(alg_name_m2)
alg_m2.set_header(["#{alg_name_m2}Alg/#{alg_name_m2}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .note(:ks_vertex_criteria, "K_S0 candidates are built from oppositely charged track pairs by a secondary vertex fit without PID; loose (unconstrained) vertex quality requirements and a decay length > 0.25 cm are imposed. The vertex-quality/decay-length selection is not expressible in the DSL and is applied via the vertex-fit criteria / downstream ROOT selection.")

sel_m2 = Selection.new
  .select_track {                      # charged-track quality cuts + multiplicity (K_S → π+π- plus K and π)
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"                        # at least two positive tracks
    nChrn ">=2"                        # at least two negative tracks
    nNet  "==0"                        # net charge zero
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"                         # at least one photon (the radiative γ)
  }
  .assign({:chrgp => :pip, :chrgn => :pim})   # K_S daughters assigned as pions, NO PID
  .secondary_vertex_fit([:pip, :pim]) {       # K_S0 from two oppositely charged tracks
    build_virtual_particle(:K_S0).by_minimizing_mass_difference   # closest to nominal K_S mass
    remove_used_particle_from_candidate_list
  }
  .pid(method: :probability) {                # remaining tracks: separate π/K/p and K/π/p
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    identify :kaon, against: [:pion, :proton]
  }
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) { # 4C fit of γ K_S K± π∓ (charge conjugate covered)
    nominal
    constrain_four_momentum
    chi2_cut 60
  }

alg_m2.with_decay_card(decay_card_kskpi).apply(sel_m2)
root_files_m2 = alg_m2.execute_on([psip_data, psip_incMC, exMC_kskpi])


# ---------------- Mode 3: χ_c2 → π+ π- π0 ----------------
alg_name_m3 = "ChiC2GammaPiPiPi0"
alg_m3 = Algorithm.new(alg_name_m3)
alg_m3.set_header(["#{alg_name_m3}Alg/#{alg_name_m3}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

sel_m3 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"                        # exactly one positive track
    nChrn "==1"                        # exactly one negative track
    nNet  "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=3"                         # at least three photons (radiative γ + π0 → γγ)
  }
  .pid(method: :probability) {         # separate π from K and p
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"                         # at least one π+
    npim ">=1"                         # at least one π-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # π0 from two photons
    invariant_mass_of(:gamma, :gamma).within(0.125, 0.145)   # |Mγγ - mπ0| < 10 MeV
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pi0]) { # 4C fit of γ π+ π- π0
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # additional π0 mass constraint
    invariant_mass_of(:pip, :pim).within(0.0, 3.0)                        # M(π+π-) in 0-3 GeV
    chi2_cut 60
  }

alg_m3.with_decay_card(decay_card_pipipi0).apply(sel_m3)
root_files_m3 = alg_m3.execute_on([psip_data, psip_incMC, exMC_pipipi0])