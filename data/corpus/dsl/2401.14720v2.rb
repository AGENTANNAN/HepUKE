# e+e- -> omega chic1,2 cross sections at sqrt(s)=4.308-4.951 GeV
# BESIII: arXiv:2401.14720v2
# omega -> pi+ pi- pi0, pi0 -> gamma gamma
# chicJ -> gamma J/psi, J/psi -> l+ l- (l=e,mu)

### Dataset ###
data_points = [
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

# Decay card (KKMC + psi(4260) for continuum production)
decay_card_omega_chic = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c1 PHSP;
  Enddecay
  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay
  Decay chi_c1
  1.000 gamma J/psi HELAMP 1 0 1 0 1 0 0 0 1 0;
  Enddecay
  Decay J/psi
  1.000 e+ e- VLL;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_omega_chic = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_omega_chic"
  config.events        = 100_000
  config.decay_card    = decay_card_omega_chic
  config.cross_section = :default
end

### Event selection ###
alg_omega_chic = Algorithm.new("OmegaChicAnalysis")
alg_omega_chic.set_header(["OmegaChicAnalysisAlg/OmegaChicAnalysis.h"])

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     "==2"
  nChrn     "==2"
  nNet      "==0"
end
.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=3"
end
.pid(method: :probability) do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 1.0
  nlp ">=1"
  nlm ">=1"
end
.remove([:lp <= :chrgp, :lm <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=1"
end
# 5C kinematic fit: 4-momentum + π0 mass constraint (π0 pre-reconstructed above)
.kinematic_fit([:pip, :pim, :lp, :lm, :gamma, :pi0]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg_omega_chic
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 5C kinematic fit; efficiency difference with/without taken as systematic uncertainty")
  .note(:omega_jpsi_mass_windows, "omega mass window [0.75, 0.81] GeV/c2 and J/psi mass window [3.08, 3.12] GeV/c2 applied at ROOT level (~3 sigma of signal MC shape)")
  .note(:chi2_5C_tight, "published tight chi2_5C < 60 cut not applied in BOSS DSL (Rule T3: loose chi2_cut 200); tight cut determined by optimization applied in ROOT")
  .note(:both_chic_states, "analysis covers both chi_c1 and chi_c2 signals with identical selection; decay card uses chi_c1 as representative")
  .note(:both_lepton_channels, "J/psi -> e+e- and mu+mu- both used; high-momentum lepton ID handles both; decay card uses e+e-")
  .note(:isr_correction, "ISR effect included in signal MC generation; radiative correction factor (1+delta) obtained iteratively from cross section measurement")
  .with_decay_card(decay_card_omega_chic)
  .apply(event_selection)

alg_omega_chic.execute_on(data_points + incMC_points + exMC_omega_chic)