# Test of CP Symmetry in the Neutral Decays of Lambda via J/psi -> Lambda anti-Lambda
# BESIII Collaboration, arXiv:2510.24333v1
# Data: (10087+/-44) x 10^6 J/psi events
# Decay: J/psi -> Lambda(-> n pi0) anti-Lambda(-> anti-p pi+) and c.c.
# n reconstructed via missing mass: P_n = P_Lambda - P_pi0

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_for_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0   PHSP;
    Enddecay

    Decay Lambda0
    1.0000 n0 pi0                PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+           HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_Lambda_antiLambda_CP_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_for_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
# This analysis uses partial reconstruction: the neutron is undetected.
# The anti-Lambda is reconstructed via secondary vertex fit, and pi0 via kalman fit.
# The neutron is inferred via partial_miss.

alg_name = "LambdaCP"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
             .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93
                 Vz 10.0
                 Vr 1.0
                 nChrp ">=1"
                 nChrn ">=1"
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:kaon, :pion]
                 nprm ">=1"
               }
               .remove([:prm <= :chrgn])
               .assign({chrgp: :pip, chrgn: :pim})
               .secondary_vertex_fit([:prm, :pip]) do
                 build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               end
               .select_photon {
                 tdc_emc_start 0
                 tdc_emc_end 14
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 nGam ">=2"
               }
               .kalman_kinematic_fit([:gamma, :gamma]) do
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 chi2_cut 100
                 npi0 ">=1"
               end
               # Partial reconstruction: miss the neutron (recID 1 from decay card: Lambda -> n pi0)
               # Lambda_bar and pi0 are pre-built; n0 is inferred from recoil
               .partial_miss([1]) do
                 best_combination_by_mass :Lambda_bar, 1.115683
                 require_recoil_mass 0.926, 0.957
               end
# Post-selection steps (ROOT):
# - Recoil mass cut on anti-Lambda (pre-selection)
# - BDT classifier for neutron suppression (BDT > 0.15) -- not expressible in DSL
# - Lambda mass region signal extraction [0.926, 0.957] GeV/c^2
# - n/anti-n separation via missing mass: P_n = P_Lambda - P_pi0
# - CP asymmetry observables from n pi0 vs anti-n pi0 decay rates
my_Algorithm.note(:bdt_classifier, "BDT > 0.15 for neutron suppression not expressible in BOSS DSL")

my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])