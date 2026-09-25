# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data      = DatasetManager.real_data.find("709_3686")    # ψ(3686) real data at 3.686 GeV
psip_incMC     = DatasetManager.inclusive_mc.find("709_3686") # ψ(3686) inclusive MC
continuum_data = DatasetManager.real_data.find("709_3650")    # 3.65 GeV continuum data for background estimation

# Decay card: ψ(3686) → γ γ J/ψ, J/ψ → e+ e-
decay_card_ee = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card: ψ(3686) → γ γ J/ψ, J/ψ → μ+ μ-
decay_card_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card: ψ(3686) → γ χ_cJ, χ_cJ → γ J/ψ  (cascade background)
decay_card_cascade = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 200k exclusive MC for J/ψ → e+ e-
exMC_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gg_jpsi_ee"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_ee
    config.cross_section   = :default
end

# 200k exclusive MC for J/ψ → μ+ μ-
exMC_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gg_jpsi_mumu"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_mumu
    config.cross_section   = :default
end

# 500k exclusive MC for the ψ(3686) → γ χ_cJ → γ γ J/ψ cascade
exMC_cascade = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gammachicJ_cascade"
    config.related_dataset = psip_data
    config.events          = 500000
    config.decay_card      = decay_card_cascade
    config.cross_section   = :default
end

### Event selection (BOSS) ###
# Both J/ψ → e+e- and J/ψ → μ+μ- modes are reconstructed from the same γγ ℓ+ℓ- final state
# and share the identical selection chain, so a single algorithm is used for both.
alg_name = "PsiPrimeGGJpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # √s = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93   # |cosθ| < 0.93
                  Vz          10.0   # |Vz| < 10 cm
                  Vr          1.0    # Vr < 1 cm in the transverse plane
                  nChrp       "==1"  # exactly one positive charged track
                  nChrn       "==1"  # exactly one negative charged track
                  nNet        "==0"  # net charge zero
                }
               # Each remaining charged track must have p > 0.8 GeV/c
               .remove(:chrgp) { condition "three_momentum_of(:chrgp) < 0.8" }
               .remove(:chrgn) { condition "three_momentum_of(:chrgn) < 0.8" }
               .select_photon {
                  tdc_emc_start     0       # 700 ns timing window start
                  tdc_emc_end       14      # 700 ns timing window end
                  angle_to_track    10.0    # ≥ 10° from the nearest charged track
                  energyThreshold_b 0.025   # E > 25 MeV in the barrel (|cosθ| < 0.8)
                  energyThreshold_e 0.050   # E > 50 MeV in the endcap (0.86 < |cosθ| < 0.92)
                  nGam              "==2"   # exactly two photon candidates
                }
               .pid(method: :probability) {
                  prob_cut 0.001   # PID probability > 0.001
                  # High-momentum tracks (p > 1.0 GeV/c) treated as leptons:
                  # electron if EMC energy > 0.6 GeV, otherwise muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"        # exactly one positive lepton
                  nlm "==1"        # exactly one negative lepton
                }
               .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
                  nominal
                  vertex_fit([2, 3])                                              # common vertex for the two leptons (indices 2,3)
                  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi) # 1C: dilepton mass → nominal J/ψ mass
                  constrain_four_momentum                                         # 4C: γγ + J/ψ → initial ψ(3686)
                  chi2_cut 200                                                    # loose BOSS χ² (tight χ²/ndof < 12 applied in ROOT)
                }

my_algorithm
  .note(:lepton_pid_criteria, "high-momentum lepton PID: a track is an electron if EMC energy > 0.6 GeV with E/p > 0.7, otherwise a muon with E/p < 0.6; the E/p criteria cannot be expressed in the DSL pid block")
  .note(:vertex_fit_quality, "kinematic-fit vertex-fit χ²/ndof < 20 required for the two lepton tracks; no DSL switch exposes the vertex-fit χ² cut")
  .with_decay_card(decay_card_ee)
  .apply(event_selection)

# Execute on real data, inclusive MC, continuum data, and the three exclusive MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC, continuum_data,
                                      exMC_ee, exMC_mumu, exMC_cascade])