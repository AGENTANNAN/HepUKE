# e+e- -> pi+ pi- omega at 24 energy points (4.0-4.6 GeV)
# 5C kinematic fit, omega -> pi+ pi- pi0, pi0 -> gamma gamma
# Helicity amplitude analysis

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_points = [
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4090"),
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
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4390"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
]

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# ConExc decay card for continuum e+e- -> pi+ pi- omega
decay_card = <<~DECAYCARD
  Decay vpho
    1 ConExc 6;
  Enddecay

  Decay vhdr
    1 pi+ pi- omega PHSP;
  Enddecay

  Decay omega
    1.0 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
    1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pipiomega"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("PipimOmegaAnalysis")
alg.set_header(["PipimOmegaAnalysisAlg/PipimOmegaAnalysis.h"])

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz  10.0
  Vr  1.0
  nChrp ">=2"     # two pi+ (one from omega, one direct)
  nChrn ">=2"     # two pi- (one from omega, one direct)
  nNet  "==0"
end
.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam ">=2"        # at least 2 photons for pi0 -> gamma gamma
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip ">=2"
  npim ">=2"
end
.remove([:pip <= :chrgp, :pim <= :chrgn])
# Reconstruct pi0 from photon pairs via Kalman fit
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=1"
end
# 5C kinematic fit: pi+ pi- pi+ pi- pi0 (4-momentum + pi0 mass constraint)
.kinematic_fit([:pip, :pip, :pim, :pim, :pi0]) do
  nominal
  constrain_four_momentum
  chi2_cut 60       # chi2_5C < 60
end
# Competing hypothesis with extra photon for background suppression
.kinematic_fit([:pip, :pip, :pim, :pim, :pi0, :gamma]) do
  constrain_four_momentum
end

alg.note(:emc_over_p_veto, "E_EMC/p < 0.9 for pion candidates from non-omega decay to suppress e+e- -> gamma omega with gamma conversion")
   .note(:Ks0_veto, "|M(pi+pi-) - M(K_S0)| veto applied: pi+ pi- invariant mass outside (0.49, 0.51) GeV/c2 for all combinations; suppresses e+e- -> K_S0 pi+ pi- pi0")
   .note(:chic0_omega_veto, "|M(pi+pi-pi0) - M(chi_c0)| veto: outside (3.39, 3.44) GeV/c2; suppresses e+e- -> chi_c0 omega")
   .note(:omega_selection, "omega candidate chosen as pi+pi-pi0 combination closest to nominal omega mass; omega signal region (0.76, 0.82) GeV/c2, sideband (0.68,0.74) and (0.84,0.90) GeV/c2")
   .note(:helicity_amplitude, "Helicity amplitude analysis with intermediate states f0(500), f0(980), f2(1270), f0(1370), b1(1235)+-, rho(1450)+-; simultaneous fit over energy points performed in ROOT")
   .note(:conexc_generator, "ConExc generator used for continuum e+e- -> pi+ pi- omega MC")
   .note(:charge_conjugate, "Charge conjugated modes implied throughout")
   .with_decay_card(decay_card)
   .apply(event_selection)

datasets = data_points + [incMC_4180, incMC_4260] + exMCs
alg.execute_on(datasets)