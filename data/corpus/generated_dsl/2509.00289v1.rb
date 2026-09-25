### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")        # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # corresponding inclusive MC sample

# Decay card for the signal process: ψ(2S) → γ χ_cJ → γ Λ anti-Λ, Λ → p π⁻, anti-Λ → anti-p π⁺
# The three radiative transitions share the same final state (γ Λ anti-Λ);
# relative weights are the PDG branching fractions of ψ(2S) → γ χ_cJ (normalised to 1).
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    0.348 gamma chi_c0 PHSP;
    0.321 gamma chi_c1 PHSP;
    0.331 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# 10M-event exclusive MC sample of the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_to_gamma_chicJ_to_LLbar"
  config.related_dataset = psip_data          # associate with the ψ(2S) real data sample
  config.events          = 10_000_000         # 10M events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GamChicJToLLbar"
gamma_chicJ_alg = Algorithm.new(alg_name)
gamma_chicJ_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
               .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy: ψ(2S) at 3.686 GeV
               .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
# --- charged track selection ---
event_selection.select_track {
                  cos_theta 0.93      # |cosθ| < 0.93
                  Vz        20.0      # |Vz| < 20 cm
                  nChrp     ">=2"     # at least 2 positively charged tracks
                  nChrn     ">=2"     # at least 2 negatively charged tracks
                  nNet      "==0"     # net charge zero
                }
               # --- particle identification (probability method) ---
               .pid(method: :probability) {
                  prob_cut 0.001                          # PID probability > 0.001
                  identify :proton, against: [:pion]      # p / anti-p identified against pions
                  nprp ">=1"                              # at least one proton
                  nprm ">=1"                              # at least one anti-proton
                }
               # --- remove the identified (anti-)protons from the generic charged-track lists ---
               .remove([:prp <= :chrgp, :prm <= :chrgn])
               # --- remaining positive/negative tracks are taken as π⁺ / π⁻ ---
               .assign({:chrgp => :pip, :chrgn => :pim})
               # --- Λ → p π⁻ secondary-vertex fit ---
               .secondary_vertex_fit([:prp, :pim]) {
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # --- anti-Λ → anti-p π⁺ secondary-vertex fit ---
               .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # --- photon selection ---
               .select_photon {
                  tdc_emc_start     0        # TDC start time
                  tdc_emc_end       14       # TDC end time
                  energyThreshold_b 0.025    # EMC barrel energy threshold  25 MeV
                  energyThreshold_e 0.050    # EMC endcap energy threshold  50 MeV
                  nGam              ">=1"    # at least one photon
                }
               # --- Λ / anti-Λ mass window [1.108, 1.123] GeV/c² ---
               .remove(:Lambda) { condition "invariant_mass_of(:Lambda) < 1.108 || invariant_mass_of(:Lambda) > 1.123" }
               .remove(:Lambda_bar) { condition "invariant_mass_of(:Lambda_bar) < 1.108 || invariant_mass_of(:Lambda_bar) > 1.123" }
               # --- kinematic fit to γ Λ anti-Λ ---
               .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) {
                  nominal                                            # nominal fit
                  constrain_four_momentum                            # 4C energy-momentum constraint
                  chi2_cut 40                                        # χ² < 40
                  invariant_mass_of(:gamma, :Lambda).out_of(1.183, 1.216)      # veto M(γΛ)   ∈ [1.183, 1.216] GeV/c²
                  invariant_mass_of(:gamma, :Lambda_bar).out_of(1.183, 1.219)  # veto M(γ anti-Λ) ∈ [1.183, 1.219] GeV/c²
                }

# BOSS-side procedure that cannot be expressed in the DSL chain
gamma_chicJ_alg
  .note(:decay_length_significance, "Λ and anti-Λ candidates are required to have a decay-length
    significance L/σ_L > 2.0, derived from the secondary-vertex fit (flight length over its error);
    the DSL secondary_vertex_fit block exposes no flight-length-significance criterion, so this
    BOSS-side V0 cut is applied on top of the fit output")

# Generate the algorithm for the process defined in the decay card and execute it
gamma_chicJ_alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = gamma_chicJ_alg.execute_on([psip_data, psip_incMC, exMC_signal])