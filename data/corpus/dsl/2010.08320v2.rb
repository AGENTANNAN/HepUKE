# e+e- -> Xi- anti-Xi+ Born cross-section measurement at 8 CMS energies (2.644-3.080 GeV)
# Single baryon tag: Xi- -> Lambda pi-, Lambda -> p pi-, anti-Xi+ inferred from recoil

decay_card_xi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Xi- anti-Xi+  PHSP;
  Enddecay

  Decay Xi-
  1.0000 Lambda pi-  PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda pi+  PHSP;
  Enddecay

  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay

  Decay anti-Lambda
  1.0000 anti-p- pi+  PHSP;
  Enddecay

  End
DECAYCARD

# Energy scan data points: use psi(4260) convention for scan data
scan_data = [
  DatasetManager.real_data.find("703_2644"),  # 2.644 GeV, 33.7 pb-1
  DatasetManager.real_data.find("703_2646"),  # 2.646 GeV, 34.0 pb-1
  DatasetManager.real_data.find("703_2900"),  # 2.900 GeV, 105 pb-1
  DatasetManager.real_data.find("703_2950"),  # 2.950 GeV, 15.9 pb-1
  DatasetManager.real_data.find("703_2981"),  # 2.981 GeV, 16.1 pb-1
  DatasetManager.real_data.find("703_3000"),  # 3.000 GeV, 15.9 pb-1
  DatasetManager.real_data.find("703_3020"),  # 3.020 GeV, 17.3 pb-1
  DatasetManager.real_data.find("703_3080"),  # 3.080 GeV, 126 pb-1
]

exMC_xi = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_xixbar"
  config.events = 100_000
  config.decay_card = decay_card_xi
  config.cross_section = :default
end

algorithm = Algorithm.new("XiXbarAnalysis")
algorithm.set_header(["XiXbarAnalysisAlg/XiXbarAnalysis.h"])
          .set_constant({"ECMS" => [:double, 3.080]})
          .note(:single_tag_method,
            "Single baryon tag: fully reconstruct Xi- -> Lambda pi- (Lambda -> p pi-); " \
            "anti-Xi+ inferred from recoiling mass M_recoil of Lambda pi- system. " \
            "Double-counting correction factor ~19% applied at ROOT level.")
          .note(:isr_correction,
            "ISR correction factor (1+delta) obtained via iteration method: " \
            "feed measured Born cross-section back into MC until convergence at 1.0% level. " \
            "Radiative correction applied in ROOT analysis.")
          .note(:upper_limit_method,
            "At threshold (2.644, 2.646 GeV) no significant signal: " \
            "90% CL upper limits via profile likelihood method.")

event_selection = Selection.new
event_selection.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=3"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  npim ">=2"
}
.remove([:prp <= :chrgp])
# Reconstruct Lambda -> p pi-
.secondary_vertex_fit([:prp, :pim]) {
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}
# Reconstruct Xi- -> Lambda pi-
.secondary_vertex_fit([:Lambda, :pim]) {
  build_virtual_particle(:Xi).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}
# Kinematic fit with Xi- and infer anti-Xi+ as recoil
.kinematic_fit([:Xi, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}
.note(:recoil_mass,
  "anti-Xi+ inferred from M_recoil = sqrt((E_cm - E_Lambda_pi)^2 - |p_cm - p_Lambda_pi|^2). " \
  "Signal yield from unbinned max-likelihood fit to M_recoil spectrum.")

algorithm.with_decay_card(decay_card_xi).apply(event_selection)
algorithm.execute_on(scan_data + exMC_xi)