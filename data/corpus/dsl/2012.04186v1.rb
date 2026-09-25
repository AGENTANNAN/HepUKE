# BESIII DSL: e+e- -> J/psi X inclusive cross sections at 69 c.m. energies
# Paper: 2012.04186v1
# Energy range: 3.645 - 3.891 GeV, KKMC generator

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ============================================================
# Major datasets at key energy points
# Use highest BOSS versions (most recent rounds)
# ============================================================
data_3650 = DatasetManager.real_data.find("709_3650")   # 3.650 GeV, 44.5 pb^-1
data_3682 = DatasetManager.real_data.find("709_3682")   # 3.682 GeV, scan
data_3686 = DatasetManager.real_data.find("709_3686")   # 3.686 GeV, 3208.5 pb^-1 (round14)
data_3773 = DatasetManager.real_data.find("712_3773")   # 3.773 GeV, 8157 pb^-1 (round16)
data_3808 = DatasetManager.real_data.find("712_3800")   # 3.800 GeV, scan (closest to 3.808)

# ============================================================
# Note: Paper uses 69 c.m. energies from 3.645 to 3.891 GeV.
# Many are from high-resolution cross-section scan data (ref 15, Chin.Phys.C 42 063001).
# The dataset table may not enumerate all 69 energy points individually.
# Additional scan energies (3.645-3.678, 3.679-3.696, 3.700-3.773, 3.776-3.891 GeV)
# must be located in the scan dataset directories and added manually.
# The 712 scan entries cover 3.800-3.885 in ~5 MeV steps.
# ============================================================

# All identifiable major datasets
major_datasets = [data_3650, data_3682, data_3686, data_3773, data_3808]

# Find corresponding inclusive MC
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")
incMC_3682 = DatasetManager.inclusive_mc.find("709_3682")
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_3808 = DatasetManager.inclusive_mc.find("712_3800")

# ============================================================
# Decay cards
# ============================================================

# J/psi -> e+ e- signal card
decay_card_jpsi_ee = <<~DECAYCARD
  Decay J/psi
    1 e+ e- PHSP;
  Enddecay
  End
DECAYCARD

# J/psi -> mu+ mu- signal card
decay_card_jpsi_mumu = <<~DECAYCARD
  Decay J/psi
    1 mu+ mu- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for J/psi -> e+ e- (at 3.686 GeV, rescaled for other energies)
sig_jpsi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_ee"
  config.related_dataset = data_3686
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi_ee
  config.cross_section   = :default
end

sig_jpsi_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_mumu"
  config.related_dataset = data_3686
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi_mumu
  config.cross_section   = :default
end

# ============================================================
# Note: Signal MC for inclusive J/psi X is generated with KKMC,
# not EvtGen. The KKMC generator simulates e+e- -> psi(3686)/psi(3770)
# with ISR effects. The decays psi -> J/psi pi pi, J/psi eta,
# J/psi pi0, gamma chi_cJ are then handled by EvtGen.
# This hybrid KKMC+EvtGen generation is configured in the BOSS
# job options, not expressible in pure EvtGen decay card syntax.
# ============================================================

# ============================================================
# Algorithm: Inclusive J/psi -> l+ l- X selection
# J/psi reconstructed via e+e- or mu+mu- pair
# ============================================================

alg_jpsi = Algorithm.new("JpsiInclusiveX", version: '00-00-01')
alg_jpsi.set_header(["JpsiInclusiveXAlg/JpsiInclusiveX.h"])
         .set_constant({ "ECMS" => [:double, 3.686] })   # nominal, overridden per energy

# ============================================================
# Note: The analysis does NOT use a conventional BOSS kinematic fit.
# Instead, the numbers of J/psi candidates are determined by fitting
# the l+l- invariant mass spectrum in ROOT after applying selection cuts.
# The selection here applies all pre-fit quality cuts.
# ============================================================

jpsi_selection = Selection.new
  .select_track do
    nChrp ">=1"
    nChrn ">=1"
    cos_theta 0.93         # |cos(theta)| < 0.93
    Vz 10.0                # |Vz| < 10 cm
    Vr 1.0                 # |Vr| < 1 cm
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025    # 25 MeV barrel
    energyThreshold_e 0.050    # 50 MeV endcap
    angle_to_track    10.0     # > 10 degrees from nearest track
    nGam              ">=1"
  end
  .pid do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.7
    nlp ">=1"; nlm ">=1"    # at least one e+ and one e- candidate
  end

# Notes for inexpressible BOSS-side procedures
alg_jpsi.note(:lepton_ep_cut, "Electron: E/p > 0.7; Muon: 0.05 < E/p < 0.35; applied in ROOT analysis after reconstruction")
alg_jpsi.note(:lepton_cos_theta, "Lepton |cos(theta)| < 0.81; angle between leptons < 179 deg; applied in ROOT")
alg_jpsi.note(:lepton_momentum, "Lepton momentum > 1 GeV and < 0.47*E_cm; applied in ROOT")
alg_jpsi.note(:pion_pid, "Pion PID: CL(pi) > CL(K) from dE/dx + TOF; applied for multi-track events via ROOT")
alg_jpsi.note(:track_mult, "Exactly 2 tracks + >=1 photon, or 3-4 charged tracks; applied in ROOT event categorization")
alg_jpsi.note(:jpsi_mass_fit, "J/psi yield from fit to l+l- invariant mass spectrum in ROOT; signal shape from MC, background = line")
alg_jpsi.note(:kkmc_generator, "Signal MC uses KKMC generator for e+e- -> psi(3686)/psi(3770) with ISR, not pure EvtGen")
alg_jpsi.note(:backgrounds, "Backgrounds: (gamma)e+e-, (gamma)mu+mu-, (gamma)tau+tau-, (gamma)DDbar, continuum, (gamma)J/psi; estimated from dedicated MC and subtracted in ROOT")
alg_jpsi.note(:efficiency, "Selection efficiency 58.8%-60.8% (energy-dependent), determined from KKMC MC with J/psi->pi+pi-, pi0pi0, eta, pi0, gamma chi_cJ decays")
alg_jpsi.note(:scan_energies, "69 c.m. energies from 3.645-3.891 GeV; 72 pb^-1 scan data (ref 15) + 44.5 pb^-1 at 3.650 + 162.8 pb^-1 at 3.686 + 2.93 fb^-1 at 3.773 + 50.5 pb^-1 at 3.808")

alg_jpsi.with_decay_card(decay_card_jpsi_ee).apply(jpsi_selection)
alg_jpsi.execute_on(major_datasets + [incMC_3650, incMC_3682, incMC_3686, incMC_3773, incMC_3808, sig_jpsi_ee, sig_jpsi_mumu])