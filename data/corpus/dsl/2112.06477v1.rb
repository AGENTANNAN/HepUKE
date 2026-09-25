DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 28 energy points from 4.085 to 4.600 GeV (Tables 1 and 2 of the paper)
scan_points = [
  DatasetManager.real_data.find("703_4090"),  # 4.0854 GeV
  DatasetManager.real_data.find("705_4130"),  # 4.1285 GeV
  DatasetManager.real_data.find("705_4160"),  # 4.1574 GeV
  DatasetManager.real_data.find("703_4180"),  # 4.1780 GeV
  DatasetManager.real_data.find("703_4190"),  # 4.1886 GeV
  DatasetManager.real_data.find("703_4200"),  # 4.1989 GeV
  DatasetManager.real_data.find("703_4210"),  # 4.2092 GeV
  DatasetManager.real_data.find("703_4220"),  # 4.2171 GeV
  DatasetManager.real_data.find("703_4230"),  # 4.2263 GeV
  DatasetManager.real_data.find("703_4237"),  # 4.2357 GeV
  DatasetManager.real_data.find("703_4246"),  # 4.2438 GeV
  DatasetManager.real_data.find("703_4260"),  # 4.2580 GeV
  DatasetManager.real_data.find("703_4270"),  # 4.2668 GeV
  DatasetManager.real_data.find("703_4280"),  # 4.2777 GeV
  DatasetManager.real_data.find("705_4290"),  # 4.2879 GeV
  DatasetManager.real_data.find("703_4310"),  # 4.3079 GeV
  DatasetManager.real_data.find("705_4315"),  # 4.3121 GeV
  DatasetManager.real_data.find("705_4340"),  # 4.3374 GeV
  DatasetManager.real_data.find("703_4360"),  # 4.3583 GeV
  DatasetManager.real_data.find("705_4380"),  # 4.3774 GeV
  DatasetManager.real_data.find("703_4390"),  # 4.3874 GeV
  DatasetManager.real_data.find("705_4400"),  # 4.3964 GeV
  DatasetManager.real_data.find("703_4420"),  # 4.4156 GeV
  DatasetManager.real_data.find("705_4440"),  # 4.4362 GeV
  DatasetManager.real_data.find("703_4470"),  # 4.4671 GeV
  DatasetManager.real_data.find("703_4530"),  # 4.5271 GeV
  DatasetManager.real_data.find("703_4575"),  # 4.5745 GeV
  DatasetManager.real_data.find("703_4600"),  # 4.5995 GeV
]

scan_incMCs = scan_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

# Decay card: e+ e- -> D*+ D*- (signal channel 1)
# D*+ -> pi+ D0, D0 -> K- pi+; D*- not reconstructed (inferred from recoil)
decay_card_dstar_dstar = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ D*- PHSP;
    Enddecay
    Decay D*+
    1.000 pi+ D0 PHSP;
    Enddecay
    Decay D*-
    1.000 pi- anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card: e+ e- -> D*+ D- (signal channel 2)
# D- not reconstructed (inferred from recoil)
decay_card_dstar_d = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ D- PHSP;
    Enddecay
    Decay D*+
    1.000 pi+ D0 PHSP;
    Enddecay
    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay
    End
DECAYCARD

exMC_dstar_dstar = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_dstar_dstar"
  config.events        = 100_000
  config.decay_card    = decay_card_dstar_dstar
  config.cross_section = :default
end

exMC_dstar_d = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_dstar_d"
  config.events        = 100_000
  config.decay_card    = decay_card_dstar_d
  config.cross_section = :default
end

# Common algorithm: reconstruct D*+ -> pi+ D0, D0 -> K- pi+; other side from recoil
alg_dstar = Algorithm.new("DstarCrossSection")
alg_dstar.set_header(["DstarCrossSectionAlg/DstarCrossSection.h"])
          .note(:pion_momentum_split,
            "Pions split by momentum: pi_L (p < 0.3 GeV/c, from D*+ -> pi+ D0) " \
            "and pi_H (p > 0.3 GeV/c, from D0 -> K- pi+). " \
            "At least one pi_L and one pi_H required per event.")
          .note(:pid_kaon,
            "Kaon candidates: P(K) > P(pi) and P(K) > 0.001, with momentum > 0.3 GeV/c.")
          .note(:pid_pion,
            "Pion candidates: P(pi) > P(K) and P(pi) > 0.001.")
          .note(:d0_mass_window,
            "D0 mass window: M(K- pi+_H) within +/- 3 sigma of known D0 mass " \
            "(1845.4 < M < 1885.2 MeV/c^2). sigma is the measured D0 mass resolution.")
          .note(:d0_mass_constraint,
            "Vertex kinematic fit constraining M(K- pi+_H) to known D0 mass (PDG). " \
            "If multiple candidates, combination with smallest vertex+kinematic fit chi2 selected.")
          .note(:missing_side,
            "D*- (or D-) is not exclusively reconstructed; inferred from energy-momentum conservation. " \
            "2D unbinned maximum-likelihood fit to RM(pi+_L D0) vs M(pi+_L D0) to extract " \
            "signal yields for D*+ D*- and D*+ D- simultaneously.")
          .note(:cross_section_calc,
            "Born cross section: sigma_B = N_sig / (L_int * (1+delta) * |1-Pi|^-2 * B1 * B2 * epsilon). " \
            "ISR correction factor (1+delta) obtained iteratively from KKMC.")
          .note(:multi_energy,
            "Cross sections measured at 28 c.m. energy points from 4.085 to 4.600 GeV. " \
            "Born cross section computed per energy point; systematic uncertainties include " \
            "luminosity (1.0%), tracking (1.0%/track), PID (1.0%/track), branching fractions, " \
            "kinematic fit, ISR correction, fit range, and signal/background shapes.")

# Reconstruct D*+ side: K- and pi+_H for D0, pi+_L for D*+
# The D0 mass constraint uses a kinematic fit; missing side not included
sel_dstar = Selection.new
  .select_track do
    cos_theta   0.93
    Vr          1.0
    Vz          10.0
    nChrp       ">=2"
    nChrn       ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkm  ">=1"
    npip ">=2"
  end
  .kinematic_fit([:km, :pip, :pip]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_dstar.with_decay_card(decay_card_dstar_dstar).apply(sel_dstar)

all_datasets = scan_points + scan_incMCs + exMC_dstar_dstar + exMC_dstar_d
alg_dstar.execute_on(all_datasets)