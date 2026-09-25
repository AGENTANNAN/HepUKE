# 2508.11276v2: e+e- -> p K- K- anti-Xi+ Born cross sections
# Partial reconstruction: Xi+ from missing mass MM(pK-K-)
# 39 energies 3.51-4.95 GeV, 20 fb^-1

DatasetManager.load_real_data("config/BES3_dataset.md")
DatasetManager.load_inclusive_mc("config/BES3_incMC.md")

# Energy scan points (3.5-4.9 GeV)
scan_points = [
  DatasetManager.real_data.find("703_3810"),
  DatasetManager.real_data.find("703_3900"),
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4090"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4245"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4390"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
]

incMC_points = scan_points.map { |d|
  DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
}

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 p+ K- K- anti-Xi+ PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pKK_Xi"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("pKKXiFinder")
alg.set_header(["pKKXiFinderAlg/pKKXiFinder.h"])
    .set_constant({ "ECMS" => [:double, 4.260] })

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion, :proton]
    nprp "==1"; nkm "==2"
  end
  .partial_miss([4]) do
    require_recoil_mass 1.28, 1.38
  end

alg.with_decay_card(decay_card).apply(event_selection)
  .note(:extra_track_veto, "Events with additional identified tracks beyond the required 1p + 2K- are discarded. Applied via nChrp '==1' and nChrn '==2' constraints.")
  .note(:signal_extraction, "Xi+ signal extracted from missing mass MM(pK-K-) distribution. Signal shape: MC-convolved Gaussian. Background: Chebyshev polynomial. Applied in ROOT analysis.")
  .note(:vertex_fit, "Vertex fit chi2 < 100 for all 3 tracks applied in ROOT, not expressible in BOSS DSL.")
  .note(:dataset_validation, "Exact list of 39 energy points 3.51-4.95 GeV to be validated against paper. Current list is representative.")

alg.execute_on(scan_points + incMC_points + sig_mc)