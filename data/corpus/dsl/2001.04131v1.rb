### Dataset preparation ###
# e+e- -> K+K-pi0pi0 PWA at 10 energy points 2.000-2.644 GeV, total 300 pb-1
# R-scan data under BOSS 713; uses ConExc generator for ISR correction

data_2000  = DatasetManager.real_data.find("713_Rscan_2000")
data_2100  = DatasetManager.real_data.find("713_Rscan_2100")
data_2125  = DatasetManager.real_data.find("713_Rscan_2125")
data_2175  = DatasetManager.real_data.find("713_Rscan_2175")
data_2200  = DatasetManager.real_data.find("713_Rscan_2200")
data_2232  = DatasetManager.real_data.find("713_Rscan_2232")
data_2309  = DatasetManager.real_data.find("713_Rscan_2309")
data_2386  = DatasetManager.real_data.find("713_Rscan_2386")
data_2396  = DatasetManager.real_data.find("713_Rscan_2396")
data_2644  = DatasetManager.real_data.find("713_Rscan_2644")

scan_points = [data_2000, data_2100, data_2125, data_2175, data_2200,
               data_2232, data_2309, data_2386, data_2396, data_2644]

# ConExc decay card: e+e- -> K+K- 2pi0 (mode 15)
# 'Particle vpho' omitted; DSL auto-injects per-point sqrt(s)
decay_card_KKpipi = <<~DECAYCARD
  Decay vpho
  1 ConExc 15;
  Enddecay
  Decay vhdr
  1 K+ K- pi0 pi0 PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for K+K-pi0pi0 (multi-energy scan)
exMC_KKpipi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_KKpipi0pi0_conexc"
  config.events        = 200_000
  config.decay_card    = decay_card_KKpipi
  config.cross_section = :default
end

### Algorithm: e+e- -> K+K-pi0pi0 (PWA input selection) ###
alg_KKpipi = Algorithm.new("KKpipi0pi0")
alg_KKpipi.set_header(["KKpipi0pi0Alg/KKpipi0pi0.h"])
           .set_constant({"ECMS" => [:double, 2.12655]})

sel_KKpipi = Selection.new
sel_KKpipi.select_track {
             cos_theta 0.93
             Vz 10.0
             Vr 1.0
             nChrp ">=1"
             nChrn ">=1"
             nNet "==0"
           }
           .select_photon {
             tdc_emc_start 0
             tdc_emc_end 14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track 10.0
             nGam ">=4"
           }
           .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion]
             nkp "==1"
             nkm "==1"
           }
           # Reconstruct two pi0 candidates from photon pairs
           .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0 ">=2"
           }
           # 6C kinematic fit: energy-momentum + two pi0 mass constraints
           .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) {
             nominal
             constrain_four_momentum
             chi2_cut 80
           }

alg_KKpipi
  .note(:pwa, "Partial-wave analysis performed with GPUPWA framework on surviving candidates; intermediate states disentangled via unbinned maximum likelihood fit using MINUIT. PWA amplitudes and Born cross section extraction are ROOT-level procedures, not expressible in BOSS DSL.")
  .note(:conexc_generator, "Signal MC generated with ConExc package incorporating higher-order ISR correction. ISR correction factor and vacuum polarization factor obtained from generator log for cross section calculation.")
  .note(:energy_groups, "Group I data (2.000-2.232 GeV) use same intermediate processes as sqrt(s)=2.125 GeV fit; Group II data (2.309-2.644 GeV) use same processes as sqrt(s)=2.396 GeV fit.")
  .note(:background, "Backgrounds from e+e- -> e+e-, mu+mu-, gamma gamma generated with Babayaga; e+e- -> hadrons with Luarlw; two-photon events with Bestwogam. Studies indicate negligible backgrounds after selection.")
  .note(:phi2170, "Structure observed at M=2126.5 MeV with width 106.9 MeV, identified as phi(2170). Cross sections measured for subprocesses: K+(1460)K-, K1+(1400)K-, K1+(1270)K-, K*+(892)K*-(892).")
  .with_decay_card(decay_card_KKpipi)
  .apply(sel_KKpipi)

# Execute on all scan points
alg_KKpipi.execute_on(scan_points)