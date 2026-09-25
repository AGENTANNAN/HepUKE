# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")      # ψ(2S)/ψ(3686) real data sample
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC

# Decay cards: ψ(2S) → γ χ_cJ, χ_cJ → p pbar K_S^0 K^- π^+ (J = 0,1,2), K_S^0 → π+π-
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 P2GC0;
    Enddecay

    Decay chi_c0
    1.000 p+ anti-p- K_S0 K- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 p+ anti-p- K_S0 K- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c2
    1.000 p+ anti-p- K_S0 K- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (1M events each), one per χ_cJ; same final state and selection
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_gammachic0_ppbarKSKpi"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card_chic0
  config.cross_section = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_gammachic1_ppbarKSKpi"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card_chic1
  config.cross_section = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_gammachic2_ppbarKSKpi"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card_chic2
  config.cross_section = :default
end

### Event selection (BOSS) ###
# The three χ_cJ share the identical final state (γ p pbar K_S0 K^- π+) and selection,
# so a single Algorithm instance serves all of them.
alg_name = "ChicJToPPbarKSKPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93    # |cosθ| < 0.93
                  Vz        10.0    # |Vz| < 10 cm
                  Vr        1.0     # Vr < 1 cm
                  nChrp     ">=3"   # at least 3 positive tracks
                  nChrn     ">=3"   # at least 3 negative tracks
                  nNet      "==0"   # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025  # 25 MeV (barrel)
                  energyThreshold_e 0.050  # 50 MeV (endcap)
                  nGam              ">=1"  # at least one photon
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # protons vs K/π
                  identify :kaon, against: [:pion]           # kaons vs π
                  nprp ">=1"
                  nprm ">=1"
                }
               .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])  # remove identified p/pbar/K
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks as π+/π−
               .secondary_vertex_fit([:pip, :pim]) {       # reconstruct K_S0 from π+π−
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .kinematic_fit([:gamma, :prp, :prm, :K_S0, :km, :pip]) {  # nominal: γ p pbar K_S0 K- π+
                  nominal
                  constrain_four_momentum
                  chi2_cut 50
                }
               .kinematic_fit([:gamma, :prp, :prm, :K_S0, :kp, :pim]) {  # charge conjugate: γ p pbar K_S0 K+ π−
                  constrain_four_momentum
                  chi2_cut 50
                }

# Generate the algorithm for the shared final state (one decay card defines the variables).
my_algorithm.with_decay_card(decay_card_chic0).apply(event_selection)
# Execute on real data, inclusive MC and the three exclusive MC samples.
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])