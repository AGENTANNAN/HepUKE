# BESIII DSL: psi(3686) -> gamma chi_cJ -> gamma Lambda anti-Lambda
# Lambda -> p+ pi-, anti-Lambda -> anti-p- pi+
# Helicity amplitude and branching fraction measurement
# arXiv:2509.00289v1

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_cJ PHSP;
    Enddecay

    Decay chi_cJ
    1.000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_chicJ_Lambda_antiLambda"
  config.related_dataset = psip_data
  config.events = 10_000_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

alg = Algorithm.new("PsipToGammaChiCJToGammaLambdaAntiLambda")
alg.set_header(["PsipToGammaChiCJToGammaLambdaAntiLambdaAlg/PsipToGammaChiCJToGammaLambdaAntiLambda.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93
                 Vz 20.0
                 nChrp ">=2"
                 nChrn ">=2"
                 nNet "==0"
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:pion]
                 nprp ">=1"
                 nprm ">=1"
               }
               .remove([:prp <= :chrgp])
               .remove([:prm <= :chrgn])
               .assign({ :chrgp => :pip, :chrgn => :pim })
               .secondary_vertex_fit([:prp, :pim]) {
                 build_virtual_particle(:Lambda).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               .secondary_vertex_fit([:prm, :pip]) {
                 build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               .select_photon {
                 tdc_emc_start 0
                 tdc_emc_end 14
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 nGam ">=1"
               }
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 40
                 invariant_mass_of(:gamma, :Lambda).out_of(1.183, 1.216)
                 invariant_mass_of(:gamma, :Lambda_bar).out_of(1.183, 1.219)
               }

alg.note(:lambda_mass_window, "Lambda and anti-Lambda invariant mass required to be in [1.108, 1.123] GeV/c^2")
   .note(:decay_length, "Lambda(anti-Lambda) decay length divided by its error required > 2.0 to suppress non-Lambda backgrounds")
   .note(:pwa_fit, "Partial wave analysis performed on selected events using TF-PWA; chi_cJ masses and widths fixed to PDG values except chi_c0 width floated")
   .with_decay_card(decay_card_signal)
   .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC_signal])