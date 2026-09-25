# 2203.02682v2: e+e- -> omega pi0 and omega eta cross sections at 34 energy points (3.773-4.701 GeV)
# Ordinary analysis, multi-energy scan, KKMC generator
# Two independent signal channels (Rule T1)

HEADERS = ['kkmpw_evtgen/EvtGenMy/KKMCParticle.h']

# Decay cards per signal mode
DECAY_CARD_OMEGA_PI0 = <<~DECAY
Decay psi(4260)
1.0 omega pi0 PHSP;
Enddecay
Decay omega
1.0 pi+ pi- pi0 PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
DECAY

DECAY_CARD_OMEGA_ETA = <<~DECAY
Decay psi(4260)
1.0 omega eta PHSP;
Enddecay
Decay omega
1.0 pi+ pi- pi0 PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
Decay eta
1.0 gamma gamma PHSP;
Enddecay
DECAY

# --- Real data: 34 energy points ---
real_datasets = DatasetManager.load_real_data.find(
  # BOSS 712
  '712_3773',
  # BOSS 703 — round06 scan
  '703_3872', '703_3900',
  # BOSS 703 — round04
  '703_4009',
  # BOSS 705
  '705_4130', '705_4160',
  # BOSS 703 — round10
  '703_4180', '703_4190', '703_4200', '703_4210', '703_4220',
  # BOSS 703 — round06
  '703_4230', '703_4237', '703_4246',
  '703_4260', '703_4270', '703_4280',
  # BOSS 705
  '705_4290', '705_4315', '705_4340',
  # BOSS 703 — round06
  '703_4360', '703_4390', '703_4420', '703_4470', '703_4530',
  # BOSS 703 — round07
  '703_4600',
  # BOSS 706
  '706_4610', '706_4620', '706_4640', '706_4660', '706_4680', '706_4700'
)

# Corresponding inclusive MC
inc_datasets = DatasetManager.load_inclusive_mc.find(
  '712_3773',
  '703_4180', '703_4190', '703_4200', '703_4210', '703_4220',
  '703_4230', '703_4237', '703_4246',
  '703_4260', '703_4270', '703_4280',
  '703_4360', '703_4420',
  '703_4600',
  '706_4620', '706_4640', '706_4660', '706_4680', '706_4700'
)

# =====================================================================
# Channel 1: e+e- -> omega pi0
# omega -> pi+ pi- pi0, pi0(omega) -> gamma gamma, pi0(bachelor) -> gamma gamma
# =====================================================================

exclusive_mc_omega_pi0 = DatasetManager.create_exclusive_mc do |c|
  c.decay_card DECAY_CARD_OMEGA_PI0
  c.related_dataset real_datasets
end

alg_omega_pi0 = Algorithm.new("ee_to_omega_pi0", version: '00-00-01')
alg_omega_pi0.set_header(HEADERS)

sel_omega_pi0 = Selection.new

# Step 1: Charged tracks
sel_omega_pi0.select_track do |t|
  t.nChrp 1
  t.nChrn 1
  t.nTot 2
  t.nNet 0
  t.cos_theta_range(-0.93, 0.93)
  t.Vz 10.0
  t.Vr 1.0
end

# Step 2: Photon selection (>= 4 photons)
sel_omega_pi0.select_photon do |p|
  p.nGam 4
  p.energy_threshold_barrel 0.025
  p.energy_threshold_endcap 0.050
  p.angle_to_track 10.0
end

# Step 3: Pion PID
sel_omega_pi0.pid(method: :probability) do |pid|
  pid.identify(:pip, :pim, against: [:kp, :km, :prp, :prm])
  pid.prob_cut 0.001
end

# Step 4: kalman_kinematic_fit for pi0 from omega -> pi+ pi- pi0 (resonance pi0)
sel_omega_pi0.kalman_kinematic_fit([:gamma, :gamma]) do |kf|
  kf.mass_constraint :pi0
end

sel_omega_pi0.assign(:pi0_resonance, from: [:pi0])

# Step 5: kalman_kinematic_fit for bachelor pi0 -> gamma gamma
sel_omega_pi0.kalman_kinematic_fit([:gamma, :gamma]) do |kf|
  kf.mass_constraint :pi0
end

sel_omega_pi0.assign(:pi0_bachelor, from: [:pi0])

# Step 6: 6C kinematic fit
# 4-momentum conservation + pi0 resonance mass + pi0 bachelor mass
sel_omega_pi0.kinematic_fit([:pip, :pim, :pi0_resonance, :pi0_bachelor]) do |kf|
  kf.invariant_mass_of(:pi0_resonance).constrain_to_nominal_mass_of(:pi0)
  kf.invariant_mass_of(:pi0_bachelor).constrain_to_nominal_mass_of(:pi0)
  kf.nominal
  kf.chi2_cut 200
end

# Step 7 (competing hypothesis): 5C kinematic fit for background suppression
# Same as 6C but without bachelor pi0 mass constraint
sel_omega_pi0.kinematic_fit([:pip, :pim, :pi0_resonance, :pi0_bachelor]) do |kf|
  kf.invariant_mass_of(:pi0_resonance).constrain_to_nominal_mass_of(:pi0)
  kf.chi2_cut 200
end

alg_omega_pi0.with_decay_card(DECAY_CARD_OMEGA_PI0).apply(sel_omega_pi0)
alg_omega_pi0.execute_on([real_datasets, inc_datasets, exclusive_mc_omega_pi0].flatten)

# =====================================================================
# Channel 2: e+e- -> omega eta
# omega -> pi+ pi- pi0, pi0 -> gamma gamma, eta -> gamma gamma
# =====================================================================

exclusive_mc_omega_eta = DatasetManager.create_exclusive_mc do |c|
  c.decay_card DECAY_CARD_OMEGA_ETA
  c.related_dataset real_datasets
end

alg_omega_eta = Algorithm.new("ee_to_omega_eta", version: '00-00-01')
alg_omega_eta.set_header(HEADERS)

sel_omega_eta = Selection.new

# Step 1: Charged tracks
sel_omega_eta.select_track do |t|
  t.nChrp 1
  t.nChrn 1
  t.nTot 2
  t.nNet 0
  t.cos_theta_range(-0.93, 0.93)
  t.Vz 10.0
  t.Vr 1.0
end

# Step 2: Photon selection (>= 4 photons)
sel_omega_eta.select_photon do |p|
  p.nGam 4
  p.energy_threshold_barrel 0.025
  p.energy_threshold_endcap 0.050
  p.angle_to_track 10.0
end

# Step 3: Pion PID
sel_omega_eta.pid(method: :probability) do |pid|
  pid.identify(:pip, :pim, against: [:kp, :km, :prp, :prm])
  pid.prob_cut 0.001
end

# Step 4: kalman_kinematic_fit for pi0 from omega -> pi+ pi- pi0
sel_omega_eta.kalman_kinematic_fit([:gamma, :gamma]) do |kf|
  kf.mass_constraint :pi0
end

sel_omega_eta.assign(:pi0_resonance, from: [:pi0])

# Step 5: kalman_kinematic_fit for eta -> gamma gamma
sel_omega_eta.kalman_kinematic_fit([:gamma, :gamma]) do |kf|
  kf.mass_constraint :pi0
end

sel_omega_eta.assign(:eta, from: [:pi0])

# Step 6: 6C kinematic fit
# 4-momentum conservation + pi0 mass + eta mass
sel_omega_eta.kinematic_fit([:pip, :pim, :pi0_resonance, :eta]) do |kf|
  kf.invariant_mass_of(:pi0_resonance).constrain_to_nominal_mass_of(:pi0)
  kf.invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
  kf.nominal
  kf.chi2_cut 200
end

# Step 7 (competing hypothesis): 5C kinematic fit for background suppression
sel_omega_eta.kinematic_fit([:pip, :pim, :pi0_resonance, :eta]) do |kf|
  kf.invariant_mass_of(:pi0_resonance).constrain_to_nominal_mass_of(:pi0)
  kf.chi2_cut 200
end

alg_omega_eta.with_decay_card(DECAY_CARD_OMEGA_ETA).apply(sel_omega_eta)
alg_omega_eta.execute_on([real_datasets, inc_datasets, exclusive_mc_omega_eta].flatten)

# NOTE: Post-kinematic-fit cuts applied in ROOT analysis (not expressible in BOSS DSL):
#   - omega mass window [0.7500, 0.8150] GeV/c^2 on M(pi+ pi- pi0)
#   - chi2_5C < 60 cut for background suppression
#   - theta_gamma_gamma < 1 radian for omega eta (ISR background removal)
#   - omega sideband subtraction
#   - Best combination selection by chi2_6C minimum
#   - Bachelor pi0 momentum > resonance pi0 momentum for pi0 pair disambiguation