# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # ψ(3686) inclusive MC
cont_data  = DatasetManager.real_data.find("712_3773")       # 3.773 GeV data used as the continuum background proxy (normalised by f_c = 1.386)

# Signal decay card: ψ(3686) → K⁻ Λ(1520) Ξ̄⁺ (+ c.c.),
# Λ(1520) → p K⁻, Ξ̄⁺ → Λ̄ π⁺, Λ̄ → p̄ π⁺  ⇒  final state p p̄ K⁻ K⁻ π⁺ π⁺
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 K- Lambda(1520)0 anti-Xi- PHSP;
    Enddecay

    Decay Lambda(1520)0
    1.000 p+ K- PHSP;
    Enddecay

    Decay anti-Xi-
    1.000 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC with the full decay chain (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_KLambda1520XiBar"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "KLambda1520XiBar"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})          # CMS energy of ψ(3686) in GeV
            .set_alias({"std::vector<double>" => "Vdouble"})
            # --- BOSS-side procedures that cannot be expressed in the DSL ---
            .note(:ip_track_requirement, "at least three of the charged tracks are required to originate from the interaction point (primary-track quality requirement on top of the standard Vz/Vr cuts)")
            .note(:pid_likelihood_sum, "the p, p̄, two K⁻ and two π⁺ hypotheses are assigned by choosing, among all possible combinations, the one that maximises the sum of the PID likelihoods")
            .note(:kaon_assignment, "K₁⁻ is defined as the kaon whose invariant mass with the proton M(pK⁻) is closest to the nominal Λ(1520) mass (1.5186 GeV/c²); K₂⁻ is the remaining kaon")
            .note(:xi_recoil_tagging, "Ξ̄⁺ is tagged through the recoil mass RM(p K₁⁻ K₂⁻) (mass against the p, K₁⁻ and K₂⁻): signal region 1.31-1.34 GeV/c², sidebands 1.242-1.272 and 1.372-1.402 GeV/c²; evaluated from the fitted four-momenta")
            .note(:lambda1520_window, "Λ(1520) selected with M(p K₁⁻) in the signal region 1.50-1.54 GeV/c², sidebands 1.43-1.47 and 1.57-1.61 GeV/c²")
            .note(:continuum_scaling, "continuum background taken from the 3.773 GeV data and scaled to the ψ(3686) data with f_c = 1.386")

event_selection = Selection.new
event_selection.select_track {              # Charged track selection
    cos_theta 0.93      # |cosθ| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm
    nChrp     "==3"     # exactly three positively charged tracks
    nChrn     "==3"     # exactly three negatively charged tracks
    nNet      "==0"     # net charge zero
  }
  .pid(method: :probability) {              # PID with probability method
    prob_cut 0.001                          # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]   # one p and one p̄ against kaons/pions
    identify :km,     against: [:pion]          # two K⁻ against pions
    identify :pip,    against: [:kaon]          # two π⁺ against kaons
    nprp "==1"
    nprm "==1"
    nkm  "==2"
    npip "==2"
  }
  .secondary_vertex_fit([:prm, :pip]) {     # reconstruct Λ̄ → p̄ π⁺
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list   # remove the Λ̄ daughters from the candidate lists
  }
  .kinematic_fit([:Lambda_bar, :prp, :km, :km, :pip]) {   # 4C fit to the six charged particles p p̄ K⁻ K⁻ π⁺ π⁺
    nominal                                 # nominal fit (4-momenta from this fit are kept)
    constrain_four_momentum                 # 4-momentum constraint against the CMS energy
    chi2_cut 200                            # loose χ² cut, tight cut applied in the ROOT analysis
  }

# Generate the algorithm for the specified signal channel
my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, the continuum sample and the exclusive signal MC
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, cont_data, exMC_signal])