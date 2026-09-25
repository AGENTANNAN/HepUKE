### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # Corresponding inclusive MC sample

# Decay card for ψ(3686) → Ξ⁻ Ξ̄⁺ (Ξ⁻ → Λ π⁻, Λ → p π⁻ and charge conjugate)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi- anti-Xi+      PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi-       PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+  PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-            HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+       HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events for the full Ξ⁻Ξ̄⁺ hyperon decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_XiXibar"
  config.related_dataset = psip_data          # anchor simulation conditions to the ψ(3686) data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "XiXibar"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})          # √s = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {              # final state p p̄ π⁻ π⁻ π⁺ π⁺
                  cos_theta 0.93            # |cosθ| < 0.93
                  nChrp     ">=3"           # at least 3 positive tracks (p, π⁺, π⁺)
                  nChrn     ">=3"           # at least 3 negative tracks (p̄, π⁻, π⁻)
                  nNet      "==0"           # net charge zero
                }
               # No Vz / Vr cut and no photon selection are applied
               .pid(method: :probability) {  # PID with the probability method
                  prob_cut 0.001             # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]  # p+ and anti-p- separated from K and π
                  nprp ">=1"                 # at least one proton
                  nprm ">=1"                 # at least one anti-proton
                }
               .remove([:prp <= :chrgp, :prm <= :chrgn])  # take protons out of the charged lists
               .assign({:chrgp => :pip, :chrgn => :pim})  # remaining tracks are treated as pions
               # Λ → p π⁻ secondary vertex, pick combination closest to the Λ nominal mass
               .secondary_vertex_fit([:prp, :pim]) {
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # Λ̄ → p̄ π⁺ secondary vertex
               .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # Ξ⁻ → Λ π⁻ secondary vertex, best combination by minimal mass difference
               .secondary_vertex_fit([:Lambda, :pim]) {
                  build_virtual_particle(:Xi_minus).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # Ξ̄⁺ → Λ̄ π⁺ secondary vertex
               .secondary_vertex_fit([:Lambda_bar, :pip]) {
                  build_virtual_particle(:Xi_plus).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # Nominal 4C kinematic fit of the Ξ⁻Ξ̄⁺ system to the CMS four-momentum
               .kinematic_fit([:Xi_minus, :Xi_plus]) {
                  nominal                  # nominal fit — corrected four-momenta are the ones stored
                  constrain_four_momentum  # 4C energy–momentum constraint to √s
                  chi2_cut 200             # χ² < 200 (loose; tightened later in ROOT)
                }

my_Algorithm
  .note(:mass_window_and_decay_length, "Λ and Ξ mass windows (±5 MeV/c² around the Λ nominal mass, ±8 MeV/c² around the Ξ nominal mass) and the decay-length > 0 requirement are applied on the secondary-vertex candidates; these criteria have no dedicated DSL construct and are handled at the ROOT level.")
  .with_decay_card(decay_card_signal).apply(event_selection)

root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])