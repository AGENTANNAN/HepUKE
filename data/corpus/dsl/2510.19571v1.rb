# Evidence of Transverse Polarization of Xi0 Hyperon in psi(3686) -> Xi0 anti-Xi0
# BESIII Collaboration, arXiv:2510.19571v1
# Data: (2.712+/-0.014) x 10^9 psi(3686) events

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_for_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi0 anti-Xi0         PHSP;
    Enddecay

    Decay Xi0
    1.0000 pi0 Lambda0          PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 pi0 anti-Lambda0     PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-               HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+          HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_Xi0_antiXi0_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_for_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Xi0Pol"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
             .set_constant({"ECMS" => [:double, 3.686]})
# Paper applies momentum-dependent PID: p > 0.5 GeV/c -> proton, else pion.
# The DSL PID does not support momentum-dependent classification directly;
# this is approximated with the standard probability-based proton identification.
my_Algorithm.note(:pid_momentum_dependent, "Paper: tracks with p > 0.5 GeV/c are protons, pions otherwise. DSL uses standard probability PID.")

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93
                 Vz 10.0
                 Vr 1.0
                 nChrp ">=2"
                 nChrn ">=2"
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:kaon, :pion]
                 nprp ">=1"
                 nprm ">=1"
               }
               .remove([:prp <= :chrgp])
               .remove([:prm <= :chrgn])
               .assign({chrgp: :pip, chrgn: :pim})
               .secondary_vertex_fit([:prp, :pim]) do
                 build_virtual_particle(:Lambda).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               end
               .secondary_vertex_fit([:prm, :pip]) do
                 build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               end
               .select_photon {
                 tdc_emc_start 0
                 tdc_emc_end 14
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 nGam ">=4"
               }
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma, :gamma]) do
                 nominal
                 constrain_four_momentum
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 chi2_cut 200
               end
# Post-kinematic-fit steps (handled in ROOT):
# - Xi0/anti-Xi0 reconstruction via mass minimization of pi0+Lambda combinations
# - |M(pi0 Lambda) - m_Xi0| < 15 MeV/c^2 and |M(pi0 anti-Lambda) - m_antiXi0| < 15 MeV/c^2
# - Lambda mass window: |M(p pi) - m_Lambda| < 5 MeV/c^2, Lambda decay length > 0
# - Sideband subtraction (B1/B2/B3 regions) for background estimation
# - Unbinned maximum likelihood fit for polarization + CP observables

my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])