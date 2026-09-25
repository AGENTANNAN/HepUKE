# e+e- -> D0 D0bar and D+ D- at sqrt(s)=3.80-4.95 GeV
# BESIII: arXiv:2402.03829v3
# Single-tag technique: reconstruct one D, infer antiparticle from recoil mass
# NOT BOSS TagAnalysis: ordinary analysis with recoil mass method (Rule T1: separate algorithms for two modes)

### Dataset — XYZ and R-scan data ###
# Multi-energy scan: typical XYZ energies from the data tables
data_points = [
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
incMC_points = data_points.map { |d|
  begin
    DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
  rescue
    nil
  end
}.compact

# Decay card Mode I: D0 D0bar, tag D0 -> K- pi+ pi+ pi-
decay_card_D0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 D0 anti-D0 VSS;
  Enddecay
  Decay D0
  1.000 K- pi+ pi+ pi- PHSP;
  Enddecay
  Decay anti-D0
  1.000 anything PHSP;
  Enddecay
  End
DECAYCARD

# Decay card Mode II: D+ D-, tag D+ -> K- pi+ pi+
decay_card_Dp = <<~DECAYCARD
  Decay psi(4260)
  1.000 D+ D- VSS;
  Enddecay
  Decay D+
  1.000 K- pi+ pi+ PHSP;
  Enddecay
  Decay D-
  1.000 anything PHSP;
  Enddecay
  End
DECAYCARD

exMC_D0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_D0D0bar"
  config.events        = 50_000
  config.decay_card    = decay_card_D0
  config.cross_section = :default
end

exMC_Dp = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_DpDm"
  config.events        = 50_000
  config.decay_card    = decay_card_Dp
  config.cross_section = :default
end

### Mode I: D0 D0bar, tag via D0 -> K- pi+ pi+ pi- (Rule T1: separate Algorithm) ###
alg_D0 = Algorithm.new("D0D0barSingleTag")
alg_D0.set_header(["D0D0barSingleTagAlg/D0D0barSingleTag.h"])

sel_D0 = Selection.new
sel_D0.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     ">=2"
  nChrn     ">=2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
  nkm  ">=1"
  npip ">=2"
  npim ">=1"
end
# Partial reconstruction: tag D0 via K-pi+pi+pi-, anti-D0 inferred from recoil mass
# Use partial_miss to miss anti-D0 (recID 2, the untagged antimeson)
.partial_miss([2]) do
  best_combination_by_mass :D0, 1.86484
  require_recoil_mass 1.80, 1.95
end

alg_D0
  .note(:single_tag_technique, "single-tag technique: tag D0 via K-pi+pi+pi- with mass window 14 MeV/c2 around nominal D0 mass (~3 sigma); anti-D0 inferred from recoil mass M_D_recoil; NOT BOSS TagAnalysis")
  .note(:recoil_mass_correction, "recoil mass corrected as M_D_recoil + M_D_tag - m_D (where M_D_tag is invariant mass of K-pi+pi+pi-, m_D is nominal D mass); correction applied at ROOT level")
  .note(:charge_conjugate, "charge-conjugate mode D0bar->K+pi-pi-pi+ included by default; factor of 2 in cross section formula accounts for this")
  .note(:signal_yield_fit, "signal yield extracted by extended maximum likelihood fit to M_D_recoil spectrum [1.80, 1.95] GeV/c2; signal shape = MC shape convolved with Gaussian; background = 2nd-order polynomial")
  .note(:isr_vss, "ISR effect included in KKMC generator with VSS model; ISR correction factor obtained through QED calculation iterated until convergence")
  .with_decay_card(decay_card_D0)
  .apply(sel_D0)

### Mode II: D+ D-, tag via D+ -> K- pi+ pi+ (Rule T1: separate Algorithm) ###
alg_Dp = Algorithm.new("DpDmSingleTag")
alg_Dp.set_header(["DpDmSingleTagAlg/DpDmSingleTag.h"])

sel_Dp = Selection.new
sel_Dp.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     ">=2"
  nChrn     ">=1"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
  nkm  ">=1"
  npip ">=2"
end
.partial_miss([2]) do
  best_combination_by_mass :D_plus, 1.86966
  require_recoil_mass 1.80, 1.95
end

alg_Dp
  .note(:single_tag_technique, "single-tag technique: tag D+ via K-pi+pi+ with mass window 16 MeV/c2 around nominal D+ mass (~3 sigma); D- inferred from recoil mass M_D_recoil")
  .note(:recoil_mass_correction, "recoil mass corrected as M_D_recoil + M_D_tag - m_D; correction applied at ROOT level")
  .note(:charge_conjugate, "charge-conjugate mode D-->K+pi-pi- included by default; factor of 2 in cross section formula")
  .note(:signal_yield_fit, "signal yield extracted by extended maximum likelihood fit to M_D_recoil spectrum [1.80, 1.95] GeV/c2")
  .note(:isr_vss, "ISR effect included in KKMC generator with VSS model")
  .with_decay_card(decay_card_Dp)
  .apply(sel_Dp)

alg_D0.execute_on(data_points + incMC_points + exMC_D0)
alg_Dp.execute_on(data_points + incMC_points + exMC_Dp)