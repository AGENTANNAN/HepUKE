# 2212.07291v2: Observation of e+e- -> omega X(3872) at sqrt(s)=4.661-4.951 GeV
# X(3872) -> pi+pi-J/psi, J/psi -> l+l- (l=e,mu), omega -> pi+pi-pi0, pi0 -> gamma gamma
# Two topology types: 6-track (all charged) and 5-track (1 pion missed)
# Per Rule T1: separate Algorithm per topology

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 9 energy points: 4.661, 4.682, 4.699, 4.740, 4.750, 4.781, 4.843, 4.918, 4.951 GeV
energy_samples = %w[709_4661 709_4682 709_4699 709_4740 709_4750 709_4781 709_4843 709_4918 709_4951]

data_points = energy_samples.map { |s| DatasetManager.real_data.find(s) }
incMC_points = energy_samples.map { |s| DatasetManager.inclusive_mc.find(s) }

# Decay card: KKMC + psi(4260) top mother for standard exclusive MC
# e+e- -> omega X(3872), X(3872) -> pi+pi-J/psi, J/psi -> l+l-
# omega -> pi+pi-pi0, pi0 -> gamma gamma
decay_card = <<~DECAY
Decay psi(4260)
  1.000 omega X_3872 PHSP;
Enddecay
Decay X_3872
  1.000 pi+ pi- J/psi VSS;
Enddecay
Decay J/psi
  1.000 e+ e- PHSP;
Enddecay
Decay omega
  1.000 pi+ pi- pi0 PHSP;
Enddecay
Decay pi0
  1.000 gamma gamma PHSP;
Enddecay
End
DECAY

exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_omega_X3872"
  config.events        = 200_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# ============================================================
# Algorithm 1: 6-track topology (all charged)
# 1C kinematic fit: recoil mass of (4pi + ll) -> pi0 mass
# ============================================================

alg_6trk = Algorithm.new("OmegaX3872SixTrk", "00-00-01")
alg_6trk.set_header(["OmegaX3872SixTrk/OmegaX3872SixTrk.h"])
alg_6trk.set_constant({ "ECMS" => [:double, 4.661] })

sel_6trk = Selection.new

# Track selection: 6 charged tracks, net charge 0
sel_6trk.select_track do
  nTot "==6"
  nNet "==0"
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
end

# Photon selection: at least 1 photon
sel_6trk.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=1"
end

# PID: leptons from J/psi identified by EMC energy
# e: EMC > 0.8 GeV, mu: EMC < 0.4 GeV, at least one muon > 3 MUC layers
# Two highest-momentum (>1 GeV/c), opposite-charge tracks = lepton pair
sel_6trk.pid(method: :probability) do
  identify :ep, against: [:pion, :kaon]
  identify :em, against: [:pion, :kaon]
  identify :mup, against: [:electron, :pion, :kaon]
  identify :mum, against: [:electron, :pion, :kaon]
  identify_high_momentum_leptons(:ep, :em, :mup, :mum)
end

# 1C kinematic fit: miss track of pi0
# Recoiling mass of all charged tracks against e+e- constrained to pi0 mass
sel_6trk.kinematic_fit([:pip, :pim, :pip, :pim, :ep, :em]) do
  nominal
  miss_track_of :pi0
  constrain_four_momentum
  chi2_cut 15
end

alg_6trk.note(:lepton_selection,
  "Two highest-momentum opposite-charge tracks (> 1 GeV/c) identified as lepton pair. " \
  "Electrons: EMC deposited energy > 0.8 GeV. Muons: EMC deposited energy < 0.4 GeV. " \
  "At least one muon must penetrate > 3 MUC layers."
)

alg_6trk.note(:kinematic_fit,
  "1C kinematic fit: recoil mass of (4 pions + lepton pair) against e+e- constrained to pi0 mass. " \
  "Chi2 < 15. Uses the 4-momentum of the non-reconstructed pi0 from the fit."
)

alg_6trk.note(:mass_windows,
  "J/psi mass window: 3.07 < M(l+l-) < 3.13 GeV/c^2. " \
  "Sidebands: [2.96, 3.05] and [3.15, 3.24] GeV/c^2. " \
  "eta veto: M(pipipi0) not in [0.52, 0.58] GeV/c^2. " \
  "psi(2S) veto: M(pipiJ/psi) not in [3.680, 3.692] GeV/c^2. " \
  "omega mass window: [0.75, 0.81] GeV/c^2. " \
  "X(3872) signal region: [3.86, 3.88] GeV/c^2. " \
  "M(pipiJ/psi) defined as M(pipill) - M(ll) + m(J/psi) for resolution improvement."
)

alg_6trk.note(:signal_extraction,
  "X(3872) signal yield from unbinned ML fit to M(pipiJ/psi). " \
  "Signal shape: MC shape convoluted with Gaussian (sigma fixed from psi(2S) control sample). " \
  "Background: linear function. " \
  "Cross section per energy point from counting method in X(3872) signal region " \
  "with sideband background subtraction."
)

alg_6trk.with_decay_card(decay_card).apply(sel_6trk)

# ============================================================
# Algorithm 2: 5-track topology (one pi missed)
# 2C kinematic fit: pi0 mass + missing pi mass constraints
# ============================================================

alg_5trk = Algorithm.new("OmegaX3872FiveTrk", "00-00-01")
alg_5trk.set_header(["OmegaX3872FiveTrk/OmegaX3872FiveTrk.h"])
alg_5trk.set_constant({ "ECMS" => [:double, 4.661] })

sel_5trk = Selection.new

# Track selection: 5 charged tracks, net charge ±1
sel_5trk.select_track do
  nTot "==5"
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
end

# Photon selection: at least 2 photons for pi0 -> gamma gamma
sel_5trk.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=2"
end

# Reconstruct pi0 from gamma gamma
sel_5trk.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=1"
end

# PID: leptons identified same as 6-track
sel_5trk.pid(method: :probability) do
  identify :ep, against: [:pion, :kaon]
  identify :em, against: [:pion, :kaon]
  identify :mup, against: [:electron, :pion, :kaon]
  identify :mum, against: [:electron, :pion, :kaon]
  identify_high_momentum_leptons(:ep, :em, :mup, :mum)
end

# 2C kinematic fit:
# - gamma gamma -> pi0 mass
# - recoil mass of (3pi + pi0 + ll) -> pi mass (missing pi)
sel_5trk.kinematic_fit([:pip, :pim, :pip, :pi0, :ep, :em]) do
  nominal
  miss_track_of :pip
  constrain_four_momentum
  chi2_cut 25
end

alg_5trk.note(:lepton_selection,
  "Same lepton selection as 6-track: two highest-momentum opposite-charge tracks (> 1 GeV/c). " \
  "Electron: EMC > 0.8 GeV. Muon: EMC < 0.4 GeV. At least one muon > 3 MUC layers."
)

alg_5trk.note(:kinematic_fit,
  "2C kinematic fit: (1) M(gamma gamma) -> pi0 mass, " \
  "(2) recoil mass of (3 pions + pi0 + lepton pair) -> pi mass. " \
  "Chi2 < 25. Multiple photon pairs: choose min chi2 combination."
)

alg_5trk.note(:mass_windows,
  "Same mass windows as 6-track: J/psi [3.07, 3.13]; eta veto [0.52, 0.58]; " \
  "psi(2S) veto [3.680, 3.692]; omega [0.75, 0.81]; X(3872) signal [3.86, 3.88] GeV/c^2."
)

alg_5trk.with_decay_card(decay_card).apply(sel_5trk)

# Execute both algorithms
alg_6trk.execute_on(data_points + incMC_points + exMCs)
alg_5trk.execute_on(data_points + incMC_points + exMCs)