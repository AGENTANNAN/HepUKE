### Dataset description ###
# 22.1 fb^-1 of e+e- collision data at c.m. energies 4.008 - 4.951 GeV (XYZ energy scan region)
data_samples = [
  DatasetManager.real_data.find("703_4009"),   # 4007.6 MeV
  DatasetManager.real_data.find("703_4090"),   # 4085.5 MeV
  DatasetManager.real_data.find("705_4130"),   # 4128.5 MeV
  DatasetManager.real_data.find("705_4160"),   # 4157.4 MeV
  DatasetManager.real_data.find("703_4180"),   # 4178.0 MeV
  DatasetManager.real_data.find("703_4190"),   # 4188.8 / 4188.6 MeV
  DatasetManager.real_data.find("703_4200"),   # 4198.9 MeV
  DatasetManager.real_data.find("703_4210"),   # 4209.2 / 4207.7 MeV
  DatasetManager.real_data.find("703_4220"),   # 4218.7 / 4217.1 MeV
  DatasetManager.real_data.find("703_4230"),   # 4226.3 MeV
  DatasetManager.real_data.find("703_4237"),   # 4235.7 MeV
  DatasetManager.real_data.find("703_4245"),   # 4241.7 MeV
  DatasetManager.real_data.find("703_4246"),   # 4243.8 MeV
  DatasetManager.real_data.find("703_4260"),   # 4258.0 MeV
  DatasetManager.real_data.find("703_4270"),   # 4266.8 MeV
  DatasetManager.real_data.find("703_4280"),   # 4277.7 MeV
  DatasetManager.real_data.find("703_4310"),   # 4307.9 MeV
  DatasetManager.real_data.find("705_4290"),   # 4287.9 MeV
  DatasetManager.real_data.find("705_4315"),   # 4312.1 MeV
  DatasetManager.real_data.find("705_4340"),   # 4337.4 MeV
  DatasetManager.real_data.find("703_4360"),   # 4358.3 MeV
  DatasetManager.real_data.find("705_4380"),   # 4377.4 MeV
  DatasetManager.real_data.find("703_4390"),   # 4387.4 MeV
  DatasetManager.real_data.find("705_4400"),   # 4396.5 MeV
  DatasetManager.real_data.find("703_4420"),   # 4415.6 MeV
  DatasetManager.real_data.find("705_4440"),   # 4436.2 MeV
  DatasetManager.real_data.find("703_4470"),   # 4467.1 MeV
  DatasetManager.real_data.find("703_4530"),   # 4527.1 MeV
  DatasetManager.real_data.find("703_4575"),   # 4574.5 MeV
  DatasetManager.real_data.find("703_4600"),   # 4599.5 MeV
  DatasetManager.real_data.find("706_4610"),   # 4611.9 MeV
  DatasetManager.real_data.find("706_4620"),   # 4628.0 MeV
  DatasetManager.real_data.find("706_4640"),   # 4640.9 MeV
  DatasetManager.real_data.find("706_4660"),   # 4661.2 MeV
  DatasetManager.real_data.find("706_4680"),   # 4681.9 MeV
  DatasetManager.real_data.find("706_4700"),   # 4698.8 MeV
  DatasetManager.real_data.find("707_4740"),   # 4739.7 MeV
  DatasetManager.real_data.find("707_4750"),   # 4750.1 MeV
  DatasetManager.real_data.find("707_4780"),   # 4780.5 MeV
  DatasetManager.real_data.find("707_4840"),   # 4843.1 MeV
  DatasetManager.real_data.find("707_4914"),   # 4918.0 MeV
  DatasetManager.real_data.find("707_4946"),   # 4950.9 MeV
]
incMC_samples = data_samples.map { |ds| DatasetManager.inclusive_mc.find(ds.name) }

# Signal decay card: e+e- -> pi0 pi0 psi(3686); psi(3686) -> pi+ pi- J/psi; J/psi -> l+ l-
# Since electron+positron form the initial state, use psi(4260) as the top mother (KKMC convention).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000  pi0  pi0  psi(2S)                       PHSP;
    Enddecay

    Decay psi(2S)
    1.000  pi+  pi-  J/psi                         PHSP;
    Enddecay

    Decay J/psi
    0.500  e+   e-                                 PHOTOS  VLL;
    0.500  mu+  mu-                                PHOTOS  VLL;
    Enddecay

    Decay pi0
    1.000  gamma gamma                             PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC generated at each energy point
exMCs = DatasetManager.create_exclusive_mc_for(data_samples) do |config, ds|
  config.sample_name     = "pi0pi0_psip_signal_#{ds.name}"
  config.related_dataset = ds
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMCs.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "pi0pi0_psip"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })       # nominal CMS energy (per-run value applied)
   .note(:lepton_pi_separation,
         "Charged tracks with p > 1.1 GeV/c are assigned as leptons; tracks with p < 0.75 GeV/c " \
         "are assigned as pions. Electron and muon are separated using E/p in the EMC: " \
         "e+/e- must satisfy E/p > 0.7; mu+/mu- must satisfy E < 0.45 GeV.")
   .note(:pi0_pairing,
         "For events with more than four photons, the combination with the smallest chi2_4C is kept. " \
         "The four selected photons are combined into two pi0 candidates by minimizing " \
         "(M(g1 g2) - M_pi0)^2 + (M(g3 g4) - M_pi0)^2, with each pi0 required to satisfy " \
         "|M(gi gj) - M_pi0| < 20 MeV/c^2.")
   .note(:helix_correction,
         "Helix-parameter track correction applied to charged tracks to improve data-MC " \
         "consistency of pull distributions before evaluating the kinematic-fit efficiency.")

event_selection = Selection.new
event_selection.select_track {                     # Charged track selection
                 cos_theta 0.93                     # |cos(theta)| < 0.93
                 Vz        10.0                     # |Vz| < 10 cm along beam
                 Vr        1.0                      # |Vxy| < 1 cm in transverse plane
                 nChrp     "==2"                    # 2 positive tracks (pi+ + l+)
                 nChrn     "==2"                    # 2 negative tracks (pi- + l-)
                 nNet      "==0"                    # net charge = 0
               }
               .select_photon {                    # Photon selection
                 tdc_emc_start   0
                 tdc_emc_end     700
                 angle_to_track  10.0
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 nGam            ">=4"              # at least 4 photons (two pi0 -> gamma gamma)
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :pion, against: [:kaon, :proton]
                 identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.1,
                                                treat_as_electron_if_energy_above: 0.7
                 nlp ">=1"; nlm ">=1"
                 npip ">=1"; npim ">=1"
               }
               # 4C kinematic fit under hypothesis e+ e- -> 4 gamma pi+ pi- l+ l-; nominal fit.
               # Photons are combined into two pi0 candidates outside the fit (best pairing) and
               # subsequently mass-constrained together with J/psi in a 7C refit (see below).
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 120                       # events with chi2_4C < 120
               }
               # 7C kinematic fit: two pi0 masses + J/psi mass constrained to their nominal values
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
                 constrain_four_momentum
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:"J/psi")
                 chi2_cut 200
               }

alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = alg.execute_on(data_samples + incMC_samples + exMCs)
