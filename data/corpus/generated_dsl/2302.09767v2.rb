# =============================================================================
# ψ(3686) → Ξ0 Ξ̄0 analysis
#   Ξ0 → π0 Λ,  Ξ̄0 → π0 Λ̄,  Λ → p π−,  Λ̄ → p̄ π+,  π0 → γγ
# BOSS part: dataset preparation + event selection up to the 6C kinematic fit
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the full Ξ0 Ξ̄0 decay chain (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi0 anti-Xi0 PHSP;
    Enddecay

    Decay Xi0
    1.0000 pi0 Lambda0 PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 pi0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 5M exclusive-MC events for the complete Ξ0 Ξ̄0 decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_Xi0Xibar"
  config.related_dataset = psip_data        # anchored to the ψ(3686) data taking
  config.events          = 5_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "Xi0Xibar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])          # algorithm header
            .set_constant({ "ECMS" => [:double, 3.686] })          # ECMS = 3.686 GeV

event_selection = Selection.new
event_selection.select_track {                 # charged-track quality selection
                  cos_theta 0.93               # |cos(theta)| < 0.93
                  Vz        10.0               # |Vz| < 10 cm
                  Vr        1.0                # Vr < 1 cm
                  nChrp     ">=2"              # at least 2 positive tracks
                  nChrn     ">=2"              # at least 2 negative tracks
                }
               .select_photon {                # photon selection
                  tdc_emc_start     0          # EMC TDC window start
                  tdc_emc_end       14         # EMC TDC window end
                  energyThreshold_b 0.025      # E_gamma > 25 MeV (barrel)
                  energyThreshold_e 0.050      # E_gamma > 50 MeV (endcap)
                  angle_to_track    10.0       # angle to nearest charged track > 10 deg
                  nGam              ">=4"      # at least 4 photons (2 x pi0 -> 4 gamma)
                }
               .pid(method: :probability) {    # probability PID, prob_cut = 0.001
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]     # p+ and p-bar
                  identify :pion,   against: [:kaon, :proton]   # pi+ and pi-
                }
               # momentum requirements from the PID step: proton p > 0.5 GeV/c, pion p < 0.5 GeV/c
               .remove(:prp) { condition "three_momentum_of(:prp) < 0.5" }   # drop soft protons
               .remove(:prm) { condition "three_momentum_of(:prm) < 0.5" }   # drop soft anti-protons
               .remove(:pip) { condition "three_momentum_of(:pip) > 0.5" }   # drop hard pions
               .remove(:pim) { condition "three_momentum_of(:pim) > 0.5" }   # drop hard anti-pions
               .secondary_vertex_fit([:prp, :pim]) {                        # Lambda -> p pi-
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list                  # do not reuse the tracks
                }
               .secondary_vertex_fit([:prm, :pip]) {                        # anti-Lambda -> p-bar pi+
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # 6C kinematic fit: 4C (energy-momentum) + 2 pi0 mass constraints on 4 gamma + Lambda + Lambda_bar
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :Lambda, :Lambda_bar]) {
                  nominal                                        # nominal fit: corrected four-momenta are kept
                  constrain_four_momentum                        # 4C constraint to the CMS four-momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # first pi0 mass constraint
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # second pi0 mass constraint
                  chi2_cut 200                                   # loose chi2 < 200 (tight cut in ROOT)
                }
                # NB: the best Xi0-Xibar pair is taken from the combination with the smallest chi2
                # (default behaviour of kinematic_fit).

my_algorithm
  .note(:helix_correction, "helix-parameter correction applied to all charged tracks before the
    6C kinematic fit; the efficiency difference between running with and without the correction
    is estimated by re-running the BOSS selection on signal MC")
  .note(:efficiency_curve, "tracking / PID efficiencies are not flat; correction factors are
    derived as a function of track momentum and polar angle and applied at the ROOT level")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])