# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # inclusive ψ(2S) MC

# ---- Decay cards (EvtGen format, EvtGen particle names) ----

# Direct χ_c0 → p p̄ K+ K−
decay_card_chi_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 P2GC0;
    Enddecay

    Decay chi_c0
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Direct χ_c1 → p p̄ K+ K−
decay_card_chi_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Direct χ_c2 → p p̄ K+ K−
decay_card_chi_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c2
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Resonant χ_c1 → p̄ K+ Λ(1520), Λ(1520) → p K−
decay_card_l1520 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 anti-p- K+ Lambda(1520) PHSP;
    Enddecay

    Decay Lambda(1520)
    1.000 p+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Resonant χ_c1 → p p̄ φ, φ → K+ K−
decay_card_phi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 p+ anti-p- phi PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC samples (200k events each) ----
exMC_chi_c0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gamma_chi_c0_ppbarKK"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_chi_c0
  config.cross_section   = :default
end

exMC_chi_c1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gamma_chi_c1_ppbarKK"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_chi_c1
  config.cross_section   = :default
end

exMC_chi_c2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gamma_chi_c2_ppbarKK"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_chi_c2
  config.cross_section   = :default
end

exMC_l1520 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gamma_chi_c1_pbarK_L1520"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_l1520
  config.cross_section   = :default
end

exMC_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gamma_chi_c1_ppbar_phi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_phi
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaChiCJToPpbarKK"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})  # √s = 3.686 GeV

# All χ_cJ modes and the two resonant χ_c1 modes share the same final state
# γ p p̄ K+ K− and the same selection chain → a single Algorithm is used
# (Rule: chi_c0, chi_c1, chi_c2 share final state and selection).
event_selection = Selection.new
  .select_track {
    cos_theta 0.93      # |cosθ| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # |Vr| < 1 cm
    nChrp     "==2"     # exactly two positive charged tracks
    nChrn     "==2"     # exactly two negative charged tracks
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    tdc_emc_start     0      # TDC start
    tdc_emc_end       14     # TDC end
    angle_to_track    10.0   # > 10° from nearest charged track
    energyThreshold_b 0.08   # EMC barrel energy > 80 MeV
    energyThreshold_e 0.08   # EMC endcap energy > 80 MeV
    nGam              ">=1"  # at least one photon candidate
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and p̄ (p against K, π)
    identify :kaon,   against: [:pion, :proton] # K+ and K− (K against π, p)
    nprp "==1"   # one proton
    nprm "==1"   # one anti-proton
    nkp  "==1"   # one K+
    nkm  "==1"   # one K−
  }
  # 4C kinematic fit to γ p p̄ K+ K− (energy-momentum conservation, nominal);
  # the photon combination with smallest χ²_4C is chosen automatically.
  .kinematic_fit([:gamma, :prp, :prm, :kp, :km]) {
    nominal
    constrain_four_momentum
    # Veto intermediate Λ(1520), Λ̄(1520), φ for the direct branching fraction
    invariant_mass_of(:prp, :km).out_of(1.45, 1.59)  # |M(pK−)  − 1.52| > 0.07
    invariant_mass_of(:prm, :kp).out_of(1.45, 1.59)  # |M(p̄K+)  − 1.52| > 0.07
    invariant_mass_of(:kp,  :km).out_of(0.99, 1.05)  # |M(K+K−) − 1.02| > 0.03
    chi2_cut 200  # loose χ² cut; tight χ² and χ_cJ mass windows applied in ROOT
  }

# Attach the shared final-state decay card and render the BOSS algorithm
my_algorithm.with_decay_card(decay_card_chi_c0).apply(event_selection)

# Run on real data, inclusive MC and all exclusive MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC,
                                      exMC_chi_c0, exMC_chi_c1, exMC_chi_c2,
                                      exMC_l1520, exMC_phi])