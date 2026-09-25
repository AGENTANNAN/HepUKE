DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 23 c.m. energy scan points from 4.0 to 4.6 GeV
data_points = [
  DatasetManager.real_data.find("703_4009"),
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
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4600")
]

incMC_points = data_points.map { |dp| DatasetManager.inclusive_mc.find(dp.sample_name) }

# Decay card: e+e- -> 2(p pbar) via psi(4260) top mother (KKMC convention)
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 p+ anti-p- p+ anti-p- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for all 23 scan points
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_2ppbar"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# Algorithm for e+e- -> 2(p pbar)
alg = Algorithm.new("E2ppbar")
alg.set_header(["E2ppbarAlg/E2ppbar.h"])

# Event selection: 4 charged tracks (2p + 2pbar), PID, 3C kinematic fit
sel = Selection.new
sel.select_track {
      nChrp "==2"
      nChrn "==2"
      nNet "==0"
      cos_theta 0.93
      Vz 10.0
      Vr 1.0
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp "==2"
      nprm "==2"
    }
   .kinematic_fit([:prp, :prm, :prp, :prm]) {
      nominal
      constrain_three_momentum
      chi2_cut 60
    }

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on(data_points + incMC_points + exMC_signal)