# Dataset preparation -- multi-energy scan points for e+e- -> Sigma0 anti-Sigma0
# The analysis covers 32 center-of-mass energies from 3.50 to 4.95 GeV.
# Sample names follow convention BOSSversion_EnergyInMeV.

scan_datasets = [
  DatasetManager.real_data.find("709_3650"),
  DatasetManager.real_data.find("712_3768"),
  DatasetManager.real_data.find("712_3773"),
  DatasetManager.real_data.find("712_3780"),
  DatasetManager.real_data.find("703_3810"),
  DatasetManager.real_data.find("703_3872"),
  DatasetManager.real_data.find("703_3900"),
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4090"),
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

# ConExc decay card: mode 3 = Sigma0 anti-Sigma0
# DSL auto-injects Particle vpho <ECMS> 0.0 per energy point.
# Sigma0 -> Lambda gamma, Lambda -> p pi- (and charge conjugate)
decay_card_signal = <<~DECAYCARD
  Decay vpho
    1 ConExc 3;
  Enddecay

  Decay vhdr
    1 Sigma0 anti-Sigma0                  PHSP;
  Enddecay

  Decay Sigma0
  1.0000 Lambda gamma                     PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000 anti-Lambda gamma                PHSP;
  Enddecay

  Decay Lambda
  1.0000 p+ pi-                           PHSP;
  Enddecay

  Decay anti-Lambda
  1.0000 anti-p- pi+                      PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC: one per energy point
exMCs = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name   = "ee_Sigma0Sigmabar0"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# Include inclusive MC for each energy point
incMCs = scan_datasets.map { |ds|
  boss = ds.boss.gsub(".", "")
  DatasetManager.inclusive_mc.find("#{boss}_#{ds.sample_name}")
}.compact

# Event selection (BOSS)
alg_name = "EESigma0Sigma0bar"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     ">=2"
                  nChrn     ">=2"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]
                  nprp ">=1"
                  nprm ">=1"
                }
               .remove([:prp <= :chrgp])
               .remove([:prm <= :chrgn])
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]
                  npip ">=1"
                  npim ">=1"
                }
               .assign({:chrgp_remaining_1 => :pip, :chrgn_remaining_1 => :pim})
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=2"
                }
               .secondary_vertex_fit([:prp, :pim]) do
                 build_virtual_particle(:Lambda).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               end
               .secondary_vertex_fit([:prm, :pip]) do
                 build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               end
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) do
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               end

algorithm
  .note(:lambda_mass_window,
        "|M(p pi-) - m_Lambda| <= 5 MeV/c^2 after secondary vertex fit; " \
        "applied at ROOT level.")
  .note(:lambda_vertex_chi2,
        "Secondary vertex fit chi2 < 500 for Lambda/Lambdabar candidates; " \
        "applied at ROOT level.")
  .note(:sigma0_pairing,
        "Among Lambda Lambda_bar gamma gamma combinations, the one minimizing " \
        "DeltaM = sqrt((M(Lambda gamma1) - m_Sigma0)^2 + (M(Lambdabar gamma2) - m_Sigma0bar)^2) " \
        "is selected. Applied at ROOT level.")
  .note(:sigma0_mass_window,
        "|M(Lambda gamma) - m_Sigma0| <= 15 MeV/c^2; applied at ROOT level.")
  .note(:ecms_constant,
        "ECMS is set per energy point via create_exclusive_mc_for; " \
        "each dataset carries its own CMS energy.")
  .note(:isr_and_vp,
        "ISR correction factor (1+delta) and vacuum polarization factor 1/|1-Pi|^2 " \
        "are obtained from ConExc generator log and QED calculation respectively. " \
        "Applied at ROOT level.")
  .note(:sideband_background,
        "Background estimated from sideband regions in M(Lambda gamma) vs M(Lambdabar gamma) plane. " \
        "Applied at ROOT level.")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
algorithm.execute_on(scan_datasets + incMCs + exMCs)