# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
psip_data       = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data at 3.686 GeV
psip_incMC      = DatasetManager.inclusive_mc.find("709_3686")   # ψ(3686) inclusive MC
continuum_data  = DatasetManager.real_data.find("709_3650")      # 3.65 GeV continuum data (for 1/s-scaled subtraction)

# Decay card for the signal process ψ(3686) → ω K_S^0 K_S^0 (EvtGen format, EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 omega K_S0 K_S0 PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500,000 exclusive MC events for the signal channel
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_omega_ksks"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipOmegaKSKS"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy = 3.686 GeV
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                                # Charged track selection
    cos_theta 0.93                               # |cos(theta)| < 0.93
    Vz        10.0                               # |Vz| < 10 cm
    Vr        1.0                                # Vr < 1 cm
    nChrp     ">=3"                              # At least 3 positively charged tracks
    nChrn     ">=3"                              # At least 3 negatively charged tracks
    nNet      "==0"                              # Net charge zero
  }
  .select_photon {                               # Photon selection
    tdc_emc_start     0                          # TDC start time
    tdc_emc_end       14                         # TDC end time
    angle_to_track    10.0                       # Min angle to nearest charged track (degrees)
    energyThreshold_b 0.025                      # Barrel energy threshold (GeV)
    energyThreshold_e 0.050                      # Endcap energy threshold (GeV)
    nGam              ">=2"                      # At least 2 photons
  }
  .pid(method: :probability) {                   # Particle identification (probability method)
    prob_cut 0.001                               # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]    # Separate π+ / π- from K and p
    npip ">=3"                                   # At least 3 π+
    npim ">=3"                                   # At least 3 π-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {      # Reconstruct π0 from a photon pair (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)   # π0 mass window (GeV/c²)
    chi2_cut 25
    npi0 ">=1"                                   # At least one π0 candidate
  }
  .secondary_vertex_fit([:pip, :pim]) {          # First K_S0 → π+π- candidate
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:pip, :pim]) {          # Second (distinct) K_S0 → π+π- candidate
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:pip, :pim, :K_S0, :K_S0, :pi0]) {   # Final 4C kinematic fit
    nominal                                      # Nominal fit — corrected four-momenta used
    constrain_four_momentum                      # 4-momentum conservation to CMS energy
    invariant_mass_of(:pip, :pim, :pi0).within(0.777, 0.807)  # ω mass window (GeV/c²)
    chi2_cut 200                                 # Loose chi2 cut (tight cut applied in ROOT)
  }

# BOSS-side procedures that cannot be expressed in the formal DSL
alg.note(:ks_mass_window, "K_S0 candidates required within the mass window (0.486, 0.510) GeV/c² (3σ) — applied after the secondary vertex fit, outside the expressible DSL surface")
   .note(:ks_decay_length, "K_S0 decay length required to be > 2σ from the interaction point (displaced vertex significance cut)")
   .note(:background_veto, "combinatorial background subtracted with a 2-dimensional K_S0 K_S0 sideband method")
   .note(:efficiency_curve, "continuum contributions subtracted using the 3.65 GeV data (43.88 pb⁻¹) scaled by 1/s; signal MC modelled with a BODY3 data-driven generator using the efficiency-corrected Dalitz plot of M(ω K_S0) vs M(ω K_S0)")

# Generate the algorithm for the signal process and execute on all datasets
alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = alg.execute_on([psip_data, psip_incMC, continuum_data, exMC_signal])