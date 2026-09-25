# BOSS Ruby DSL for arXiv:2007.12872v1
# "Measurement of cross sections for e+e- -> mu+mu- at center-of-mass energies
#  from 3.80 to 4.60 GeV"
#
# Analysis: e+e- -> mu+mu- at 133 c.m. energy points from 3.8-4.6 GeV
#   Measurement of observed cross sections and Born cross sections.
#   Fits with coherent sum of continuum + 9 vector resonances (rho, omega,
#   phi, J/psi, psi(3686), psi(3770), psi(4040), psi(4160), psi(4415))
#   plus a new structure S(4220).
#   Muon ID via E/p ratio (0.05 < E/p < 0.40); Bhabha/pi+pi-/K+K-/ppbar rejection.
#   BABAYAGA generator for MC (not in standard BOSS; DSL approximates with KKMC).
#
# Key results: muonic widths and phases of psi(4040), psi(4160), psi(4415);
#   evidence for S(4220) structure at 3.9 sigma significance.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Representative subset of the 133 energy points (full scan: 133 points, 3.80-4.60 GeV).
# Listed here are the major BOSS 703/705/706/707 samples; many intermediate low-luminosity
# scan points are omitted for brevity.
scan_samples = [
  "703_3810", "703_3872", "703_3900", "703_4009", "703_4090",
  "703_4180", "703_4190", "703_4200", "703_4210", "703_4220",
  "703_4230", "703_4237", "703_4245", "703_4246", "703_4260",
  "703_4270", "703_4280", "703_4310", "703_4360", "703_4390",
  "703_4420", "703_4470", "703_4530", "703_4575", "703_4600",
  "705_4130", "705_4160", "705_4290", "705_4315", "705_4340",
  "705_4380", "705_4400", "705_4440",
  "706_4610", "706_4620", "706_4640", "706_4660", "706_4680", "706_4700",
  "707_4740", "707_4750", "707_4780", "707_4840", "707_4914", "707_4946"
]

scan_datasets = scan_samples.map { |s| DatasetManager.real_data.find(s) }
scan_incMCs   = scan_samples.map { |s| DatasetManager.inclusive_mc.find(s) }

# ===========================================================================
# Decay card: e+e- -> mu+mu- via KKMC + psi(4260)
#   BABAYAGA (the actual generator used in the paper for precise ISR/FSR)
#   is not available in the standard BOSS toolchain. The DSL approximates
#   with KKMC + VLL for the mu+mu- final state.
# ===========================================================================
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000 mu+ mu- VLL;
    Enddecay
End
DECAYCARD

# Exclusive MC for e+e- -> mu+mu- across the energy scan.
# create_exclusive_mc_for auto-injects the per-point CMS energy
# and suffixes sample names with the dataset identifier.
exMC_mumu_scan = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name = "mumu_scan_signal_mc"
  config.events = 200_000
  config.decay_card = decay_card_mumu
  config.cross_section = :default
end

# ===========================================================================
# Algorithm: e+e- -> mu+mu- selection
#   Two opposite-charge tracks with muon PID.
#   Key cuts from paper:
#     - |cos(theta)| < 0.8, Vr < 1.0 cm, |Vz| < 10.0 cm
#     - |p+| + |p-| > 0.9 * sqrt(s)   (ROOT-level; energy-dependent)
#     - 0.05 < E/p < 0.40              (muon PID; approximated via BOSS PID)
#     - 4C kinematic fit with chi2_4C < 60
#   Final signal extraction via fit to E(mu+mu-)/sqrt(s) distribution (ROOT).
# ===========================================================================
alg_mumu = Algorithm.new("EEtoMuMu")
alg_mumu.set_header(["EEtoMuMuAlg/EEtoMuMu.h"])
        .set_constant({"ECMS" => [:double, 4.200]})  # nominal mid-scan energy
        .with_decay_card(decay_card_mumu)
        .note(:full_scan_133_points, "Paper uses 133 energy points from 3.80 to 4.60 GeV. DatasetManager lists 45 major samples; remaining low-luminosity points not in dataset table")
        .note(:babayaga_not_available, "Paper uses BABAYAGA generator for e+e- -> mu+mu- MC with precise ISR/FSR/vacuum-polarization. BABAYAGA is not in standard BOSS; DSL approximates with KKMC + VLL")
        .note(:muon_pid_ep_ratio, "Muon ID uses E/p ratio: 0.05 < E(EMC)/p(track) < 0.40. Rejects Bhabha (E/p ~ 1) and pi+pi- (E/p ~ 0). DSL approximates with standard BOSS muon PID probability")
        .note(:momentum_sum_cut, "|p+| + |p-| > 0.9 * sqrt(s) cut is energy-dependent and applied per-point in ROOT after track fitting; not expressible in DSL for_each due to ECMS dependency")
        .note(:ep_ratio_energy_dependent, "E/p ratio cut efficiency is energy-dependent: pi+pi- rejection 41.5% at 3.8 GeV to 37.5% at 4.6 GeV. Handled in ROOT")
        .note(:signal_extraction, "Signal yield from double-Gaussian fit to E(mu+mu-)/sqrt(s) distribution. Peaking background (ee, pi+pi-, K+K-) subtracted via MC. ROOT-level only")
        .note(:cross_section_fit, "Observed cross sections fitted with coherent sum of continuum + 9 vector resonances + S(4220) BW. ISR correction via sampling function F(x,s). 8 solutions. ROOT-level only")
        .note(:bcdr_cross_section, "Born-continuum-dressed-resonance (BCDR) cross section derived by dividing observed xs by ISR correction factor. ROOT-level")
        .note(:offset_method_systematics, "Systematic uncertainties from correlated sources handled via offset method: cross-section values shifted simultaneously by uncertainty size. ROOT-level")

sel_mumu = Selection.new
# Charged track selection: exactly 2 tracks, opposite charge
sel_mumu.select_track do
  cos_theta 0.8       # |cos(theta)| < 0.8  (paper uses 0.8)
  Vz 10.0             # |Vz| < 10.0 cm
  Vr 1.0              # Vr < 1.0 cm
  nChrp "==1"
  nChrn "==1"
  nNet "==0"
end
# Muon PID: BOSS muon identification approximates the paper's E/p ratio cut
# (0.05 < E/p < 0.40 which rejects electrons and pions)
.pid(method: :probability) do
  prob_cut 0.001
  identify :mup, against: [:electron, :pion, :kaon, :proton]
  identify :mum, against: [:electron, :pion, :kaon, :proton]
  nmup "==1"
  nmum "==1"
end
# Four-constraint kinematic fit: e+e- -> mu+mu- hypothesis
# Paper uses chi2_4C < 60 (tighter than the standard loose 200)
.kinematic_fit([:mup, :mum]) do
  nominal
  constrain_four_momentum
  chi2_cut 60
end

alg_mumu.apply(sel_mumu)

# Execute on a representative data point + inclusive MC + exclusive MC.
# Full analysis runs over all 133 points in ROOT.
representative_data = scan_datasets.first   # 703_3810
representative_incMC = scan_incMCs.first
representative_exMC = exMC_mumu_scan.first  # first energy point from create_exclusive_mc_for

alg_mumu.execute_on([representative_data, representative_incMC, representative_exMC])