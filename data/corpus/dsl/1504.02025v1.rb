### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

### Decay cards ###
decay_card_KLambdaXi = <<~DECAYCARD
  Decay psi(2S)
  1.0000  K-  Lambda0  anti-Xi+                   PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000  anti-Lambda0  pi+                        PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                                  HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  anti-p-  pi+                             HypWK;
  Enddecay

  End
DECAYCARD

decay_card_gammaKLambdaXi = <<~DECAYCARD
  Decay psi(2S)
  1.0000  gamma  K-  Lambda0  anti-Xi+             PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000  anti-Lambda0  pi+                        PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                                  HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  anti-p-  pi+                             HypWK;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples ###
exMC_KLambdaXi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_to_K_Lambda_Xibar"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_KLambdaXi
  config.cross_section  = :default
end
exMC_KLambdaXi.save_to_config(format: :yaml, file_path: 'exMC_KLambdaXi')

exMC_gammaKLambdaXi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_to_gamma_K_Lambda_Xibar"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_gammaKLambdaXi
  config.cross_section  = :default
end
exMC_gammaKLambdaXi.save_to_config(format: :yaml, file_path: 'exMC_gammaKLambdaXi')

### Event selection: psi(3686) -> K- Lambda Xibar+ ###
alg1_name = "KLambdaXi"
alg1 = Algorithm.new(alg1_name)
alg1.set_header(["#{alg1_name}Alg/#{alg1_name}.h"])
    .set_constant({"ECMS" => [:double, 3.686]})

sel1 = Selection.new
sel1.select_track {
       cos_theta 0.93
       Vz        100.0
       Vr        10.0
       nChrp     ">=3"
       nChrn     ">=3"
       nNet      "==0"
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :pion,   against: [:kaon, :proton]
       identify :kaon,   against: [:pion, :proton]
       identify :proton, against: [:kaon, :pion]
       nprp ">=1"
       nprm ">=1"
     }
     .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .secondary_vertex_fit([:prm, :pip]) {
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .secondary_vertex_fit([:Lambda_bar, :pip]) {
       build_virtual_particle(:anti_Xi).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:km, :Lambda, :anti_Xi]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
       invariant_mass_of(:prp, :pim).within(1.110, 1.121)
       invariant_mass_of(:prm, :pip).within(1.110, 1.121)
     }

alg1.note(:kaon_ip_requirement,
          "The identified K- is required to originate from the IP: |Vxy|<1 cm and |Vz|<10 cm.")
    .note(:kaon_uniqueness,
          "Only the kaon with the highest PID CL is kept; extra kaons are treated as pions. "\
          "Same for proton/anti-proton.")
    .with_decay_card(decay_card_KLambdaXi)
    .apply(sel1)

alg1.execute_on([psip_data, psip_incMC, exMC_KLambdaXi])

### Event selection: psi(3686) -> gamma K- Lambda Xibar+ ###
alg2_name = "GammaKLambdaXi"
alg2 = Algorithm.new(alg2_name)
alg2.set_header(["#{alg2_name}Alg/#{alg2_name}.h"])
    .set_constant({"ECMS" => [:double, 3.686]})

sel2 = Selection.new
sel2.select_track {
       cos_theta 0.93
       Vz        100.0
       Vr        10.0
       nChrp     ">=3"
       nChrn     ">=3"
       nNet      "==0"
     }
     .select_photon {
       tdc_emc_start    0
       tdc_emc_end      14
       angle_to_track   10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam             ">=1"
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :pion,   against: [:kaon, :proton]
       identify :kaon,   against: [:pion, :proton]
       identify :proton, against: [:kaon, :pion]
       nprp ">=1"
       nprm ">=1"
     }
     .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .secondary_vertex_fit([:prm, :pip]) {
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .secondary_vertex_fit([:Lambda_bar, :pip]) {
       build_virtual_particle(:anti_Xi).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:gamma, :km, :Lambda, :anti_Xi]) {
       nominal
       constrain_four_momentum
       chi2_cut 100
       invariant_mass_of(:prp, :pim).within(1.110, 1.121)
       invariant_mass_of(:prm, :pip).within(1.110, 1.121)
       invariant_mass_of(:Lambda_bar, :pip).within(1.315, 1.330)
     }

alg2.note(:kaon_ip_requirement,
          "K- required to originate from the IP: |Vxy|<1 cm, |Vz|<10 cm.")
    .note(:best_photon_selection,
          "If multiple good photons, the one with the smallest 4C chi2 is selected.")
    .note(:sigma0_background_veto,
          "For chi_cJ -> K- Lambda Xibar+, additional M(gamma Lambda) > 1.21 GeV/c^2 to suppress "\
          "psi(3686) -> K- Sigma0 Xibar+; applied at ROOT stage on 4C-corrected quantities.")
    .with_decay_card(decay_card_gammaKLambdaXi)
    .apply(sel2)

alg2.execute_on([psip_data, psip_incMC, exMC_gammaKLambdaXi])
