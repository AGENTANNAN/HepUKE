# DSL for 2003.12220v2: D+ -> eta mu+ nu_mu BF measurement
# Double-tag method at psi(3770): ST D- from 6 hadronic modes,
# signal D+ -> eta mu+ nu_mu, eta -> gamma gamma
# Datasets: psi(3770), 2.93 fb^-1

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for psi(3770) -> D+ D-
decay_card_psipp = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay
  End
DECAYCARD

exMC_psipp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_DDbar_3773"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psipp
  config.cross_section   = :default
end

# TagAnalysis: ST D- (hadronic) + signal D+ -> eta mu+ nu_mu (semileptonic)
alg = TagAnalysis.new("DpEtaMuNu")
alg.set_header(["DpEtaMuNuAlg/DpEtaMuNu.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card_psipp)

# Tag side: ST D- via 6 hadronic modes
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: eta -> gamma gamma + mu+ + missing nu_mu
alg.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.missing :nu_mu
end

# Fit: 4C + eta resonance constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg.note(:st_selection, "ST D- candidates selected by deltaE and M_BC windows; deltaE in (-0.025,0.025) GeV for non-pi0 tags and (-0.055,0.045) GeV for tags with pi0. M_BC in (1.863,1.877) GeV/c^2. Best candidate per mode per charge chosen by smallest |deltaE|.")
    .note(:extra_photon_veto, "Maximum extra photon energy < 0.30 GeV and N_extra_pi0 == 0 to suppress D+ -> eta pi+ pi0 background.")
    .note(:eta_mu_inv_mass_cut, "M_eta_mu+ < 1.74 GeV/c^2 to suppress D+ -> eta pi+ background from pion-muon misidentification.")
    .note(:muon_emc_energy, "Muon candidate EMC deposited energy required in (0.105, 0.275) GeV to reduce hadron-muon misidentification.")
    .note(:umiss_signal_extraction, "U_miss = E_miss - |p_miss| used for signal extraction; D+ momentum constrained by tag side direction. Signal yield from unbinned fit to U_miss distribution.")
    .note(:q2_binning, "Partial decay rates measured in 5 q^2 intervals for form factor determination; efficiency matrix correction applied at ROOT level.")
    .note(:lepton_pid_thresholds, "Muon PID uses fixed v1 thresholds: CL_mu > 0.001, CL_mu > CL_e, CL_mu > CL_K. Signal-side lepton PID thresholds are not DSL-tunable.")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_psipp])