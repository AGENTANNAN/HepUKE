# 2202.12748v1: e+e- -> eta pi+ pi- cross section at 28 energy points (3.872-4.700 GeV)
# Ordinary analysis, multi-energy scan, KKMC generator

HEADERS = ['kkmpw_evtgen/EvtGenMy/KKMCParticle.h']

DECAY_CARD = <<~DECAY
Decay psi(4260)
1.0 eta pi+ pi- PHSP;
Enddecay
Decay eta
1.0 gamma gamma PHSP;
Enddecay
DECAY

# --- Real data: 28 energy points from 3.872 to 4.700 GeV ---
real_datasets = DatasetManager.load_real_data.find(
  '703_3872', '703_3900', '703_4009', '703_4090',
  '703_4180', '703_4190', '703_4200', '703_4210',
  '703_4220', '703_4230', '703_4237', '703_4246',
  '703_4260', '703_4270', '703_4280',
  '705_4290', '705_4315', '705_4340', '705_4380', '705_4400', '705_4440',
  '703_4360', '703_4420', '703_4470', '703_4530',
  '703_4600',
  '706_4610', '706_4620', '706_4640', '706_4660', '706_4680', '706_4700'
)

# --- Inclusive MC ---
inc_datasets = DatasetManager.load_inclusive_mc.find(
  '703_3872', '703_3900', '703_4009',
  '703_4180', '703_4190', '703_4200', '703_4210',
  '703_4220', '703_4230', '703_4237', '703_4246',
  '703_4260', '703_4270', '703_4280',
  '705_4290', '705_4315', '705_4340', '705_4380', '705_4400', '705_4440',
  '703_4360', '703_4420', '703_4470', '703_4530',
  '703_4600',
  '706_4610', '706_4620', '706_4640', '706_4660', '706_4680', '706_4700'
)

# --- Exclusive MC (multi-energy) ---
exclusive_mc = DatasetManager.create_exclusive_mc do |c|
  c.decay_card DECAY_CARD
  c.related_dataset real_datasets
end

# --- Selection ---
# NOTE: Multi-energy scan — no ECMS constant (set_constant omitted per Rule M1)

algorithm = Algorithm.new("ee_to_eta_pipi", version: '00-00-01')
algorithm.set_header(HEADERS)

event_selection = Selection.new

# Step 1: Charged tracks
event_selection.select_track do |t|
  t.nChrp 1
  t.nChrn 1
  t.nTot 2
  t.cos_theta_range(-0.93, 0.93)
  t.Vz 10.0
  t.Vr 1.0
end

# Step 2: Photon selection
event_selection.select_photon do |p|
  p.nGam 2
  p.energy_threshold_barrel 0.025
  p.energy_threshold_endcap 0.050
  p.angle_to_track 10.0
end

# Step 3: Pion PID
event_selection.pid(method: :probability) do |pid|
  pid.identify(:pip, :pim, against: [:kp, :km, :prp, :prm])
  pid.prob_cut 0.001
end

# Step 4: kalman_kinematic_fit for eta -> gamma gamma
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do |kf|
  kf.mass_constraint :pi0
end

# Assign the eta candidate
event_selection.assign(:eta, from: [:pi0])

# Step 5: 4C kinematic fit (loose chi2 for ROOT optimization)
event_selection.kinematic_fit([:eta, :pip, :pim]) do |kf|
  kf.nominal
  kf.chi2_cut 200
end

algorithm.with_decay_card(DECAY_CARD).apply(event_selection)
algorithm.execute_on([real_datasets, inc_datasets, exclusive_mc].flatten)

# NOTE: The paper applies additional background suppression cuts at the ROOT
# analysis stage (post-kinematic-fit):
#   - J/psi mass window veto on M(pi+ pi-)
#   - E/cp ratio cut for electron rejection
#   - Muon counter depth cut for muon rejection
#   - R ratio cut for continuum background suppression
# These are not expressible in the BOSS DSL and belong in the ROOT analysis.