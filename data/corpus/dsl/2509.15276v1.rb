# BESIII DSL: psi(3686) -> Lambda anti-Lambda
# Lambda -> p+ pi-, anti-Lambda -> anti-p- pi+
# Lambda transverse polarization measurement
# arXiv:2509.15276v1

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
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
  config.sample_name = "psip_Lambda_antiLambda_polarization"
  config.related_dataset = psip_data
  config.events = 5_000_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

alg = Algorithm.new("PsipToLambdaAntiLambdaPolarization")
alg.set_header(["PsipToLambdaAntiLambdaPolarizationAlg/PsipToLambdaAntiLambdaPolarization.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93
                 nChrg ">=4"
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
               .kinematic_fit([:prp, :pim, :prm, :pip]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 100
               }

alg.note(:momentum_based_pid, "Momentum-dependent PID: p < 0.6 GeV/c assigned as pion, p > 0.8 GeV/c assigned as proton; region [0.6, 0.8] excluded")
   .note(:decay_length_significance, "Lambda and anti-Lambda decay length significance L/sigma_L > 2 required")
   .note(:lambda_mass_window, "Invariant mass windows [1.108, 1.123] GeV/c^2 applied for Lambda and anti-Lambda in ROOT analysis")
   .note(:best_pair_selection, "Best Lambda-antiLambda pair chosen by minimizing (M_pi_p - m_Lambda)^2 + (M_anti-p_pi - m_Lambda)^2")
   .note(:transverse_polarization, "Lambda transverse polarization extracted via angular analysis of proton/anti-proton decay distributions")
   .note(:angular_analysis, "Angles theta_p, phi_p defined in Lambda helicity frame; asymmetry parameters measured from cos_theta distributions")
   .note(:sideband, "Lambda mass sidebands used for background subtraction in angular distributions")
   .note(:dataset_size, "448.1M psi(3686) events used for this analysis; signal yield extracted from 2D fit to M(p pi-) vs M(anti-p pi+)")
   .with_decay_card(decay_card_signal)
   .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC_signal])