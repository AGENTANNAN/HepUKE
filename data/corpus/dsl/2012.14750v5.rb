DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 16 energy scan points for CMS energy measurement (2017XYZ and 2019XYZ)
data_points = [
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("705_4440")
]

incMC_points = data_points.map { |dp| DatasetManager.inclusive_mc.find(dp.sample_name) }

# Decay card: e+e- → mu+mu- via KKMC + psi(4260) convention
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 mu+ mu- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for di-muon signal at each energy point (generated with BABAYAGA3.5 equivalent)
exMC_dimuon = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_dimuon"
  config.events        = 1_000_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("DiMuonEcmsMeasure")
alg.set_header(["DiMuonEcmsMeasureAlg/DiMuonEcmsMeasure.h"])
   .note(:cm_energy_measurement, "This analysis measures CMS energies via the di-muon process.
     M(mu+mu-) peak position is fitted with a Gaussian in (-1sigma,+1.5sigma).
     ISR/FSR mass shift correction and momentum calibration correction are applied.
     Calibration uses ISR J/psi → mu+mu- in the same data samples.")

sel = Selection.new
sel.select_track {
      nChrp "==1"
      nChrn "==1"
      nNet "==0"
      cos_theta 0.8
      Vz 10.0
      Vr 1.0
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6
      identify :muon, against: [:pion]
      nmup "==1"
      nmum "==1"
    }
   .kinematic_fit([:mup, :mum]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on(data_points + incMC_points + exMC_dimuon)