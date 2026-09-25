DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 23 energy scan points from 4.178 to 4.600 GeV
data_points = [
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
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600")
]

incMC_points = data_points.map { |dp| DatasetManager.inclusive_mc.find(dp.sample_name) }

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 pi0 X(3872) gamma PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pi0X3872gamma"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("Pi0X3872Gamma")
alg.set_header(["Pi0X3872GammaAlg/Pi0X3872Gamma.h"])
   .note(:helix_correction, "Helix parameter correction applied to charged tracks before kinematic fit")

sel = Selection.new
sel.select_track {
      nChrp ">=2"
      nChrn ">=2"
      nNet "==0"
      cos_theta 0.93
      Vz 10.0
      Vr 1.0
    }
   .select_photon {
      tdc_emc_start 0
      tdc_emc_end 700
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track 10.0
      nGam ">=3"
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6
      nlp "==1"
      nlm "==1"
    }
   .assign({:chrgp => :pip, :chrgn => :pim})
   .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
    }
   .kinematic_fit([:lp, :lm, :pip, :pim, :pi0, :gamma]) {
      nominal
      constrain_four_momentum
      chi2_cut 60
    }

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on(data_points + incMC_points + exMC_signal)