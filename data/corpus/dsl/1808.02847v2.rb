#==============================================================================
# Paper: 1808.02847v2
# Title: Evidence of a resonant structure in the e+e- -> pi+ D0 D*- cross
#        section between 4.05 and 4.60 GeV
# Journal: Phys. Rev. Lett. 122, 102002 (2019)
# Data: 84 energy points (4.05-4.60 GeV); 5 main points with >500 pb^-1
# Method: Partial reconstruction — D0 -> K-pi+, bachelor pi+;
#         D*- inferred from recoil mass RM(D0 pi+)
#==============================================================================

# --- Decay card ---
decay_card = <<~DECAY
Decay vpho
1.0 pi+ D0 D*- PHSP;
Enddecay
DECAY

# --- Algorithm ---
alg = Algorithm.new("PiD0Dstar", "00-00-01")

alg.set_header(["PiD0Dstar/PiD0Dstar.h"])
    .set_constant({ "ECMS" => [:double, 4.260] })
    .with_decay_card(decay_card)

# --- Selection ---
event_selection = Selection.new

# Charged track selection: 4 tracks (K, pi from D0, bachelor pi+, K charge opposite to pi)
event_selection.select_track do
  nTot ">=4"
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
end

# PID for kaon and pion from D0 decay
event_selection.pid(method: :probability) do
  identify :km, :pip, against: [:pip, :kp, :km, :prp, :prm, :ep, :em, :mup, :mum]
  prob_cut 0.001
end

# D0 -> K- pi+ invariant mass window and best candidate selection
# The D0 candidate is selected as the K-pi+ pair with M(K-pi+) closest to m(D0).
# After D0 selection, the bachelor pi+ is the extra charged track
# D*- inferred from recoil mass: RM_cor(D0 pi+) = RM(D0 pi+) + M(K-pi+) - m(D0)

# 5C kinematic fit constraining D0 mass and 4-momentum
# However the D*- is not directly reconstructed - partial_rec used
# Note: D*D* background vetoed by M(D0 pi+) < 2.03 GeV/c^2

# Build D0 from selected K- pi+ pair using secondary vertex
event_selection.build_virtual_particle(:D0, from: [:km, :pip])

# Partial reconstruction: reconstruct pi+ (bachelor) + D0, miss D*-
# The D*- is inferred from the recoil mass against the D0 pi+ system
# RecID mapping: 0=vpho(top,skip), 1=pi+(bachelor), 2=D0, 3=D*-(miss)
event_selection.partial_miss([3]) do
  require_recoil_mass 1.85, 2.15  # D*- mass region
end

alg.note(:d0_selection,
  "D0 candidate selected from K-pi+ pair with M(K-pi+) closest to " \
  "nominal m(D0) = 1.8648 GeV/c^2. D0 mass window: " \
  "|M(K-pi+) - m(D0)| < 15 MeV/c^2. " \
  "Only ~0.3% events have multiple D0 candidates.")

alg.note(:bachelor_pion,
  "Bachelor pi+ is the extra charged track not used in D0 reconstruction, " \
  "with charge opposite to K- from D0. If multiple candidates, choose the " \
  "one with RM_cor(D0 pi+) closest to nominal m(D*-) = 2.0103 GeV/c^2.")

alg.note(:dstar_recoil,
  "D*- not directly reconstructed. Corrected recoil mass: " \
  "RM_cor(D0 pi+) = RM(D0 pi+) + M(K-pi+) - m(D0) used. " \
  "D*- signal region: |RM_cor - DeltaM - m(D*-)| < 20 MeV/c^2. " \
  "Sideband: 1.91 < RM_cor < 1.95 GeV/c^2.")

alg.note(:dstardstar_veto,
  "e+e- -> D*D* background rejected by vetoing events with " \
  "M(D0 pi+) < 2.03 GeV/c^2.")

alg.note(:signal_fit,
  "Signal yield extracted via unbinned ML fit to RM_cor(D0 pi+) in ROOT stage. " \
  "Signal shape: MC shape convolved with Gaussian. " \
  "Background: PHSP MC (isospin partner pi+D-D*0) + linear polynomial. " \
  "Simultaneous fit to all energy points.")

alg.note(:born_cross_section,
  "Born cross section = N_obs / (L * (1+delta) * epsilon * B(D0->K-pi+)). " \
  "Resonance parameters from coherent sum of phase-space + two BW functions. " \
  "R1: M=4228.6+/-4.1+/-6.3 MeV/c^2, Gamma=77.0+/-6.8+/-6.3 MeV. " \
  "R2: M=4404.7+/-7.4 MeV/c^2, Gamma=191.9+/-13.0 MeV. " \
  "All fit work done in ROOT stage.")

alg.note(:energy_points,
  "84 energy points from 4.05 to 4.60 GeV. 5 main points with high " \
  "luminosity (>500 pb^-1): 4.2263, 4.2580, 4.3583, 4.4156, 4.5995 GeV. " \
  "79 low-luminosity points (<200 pb^-1 each). " \
  "Representative datasets: 703_4230 (1056.4 pb^-1), 703_4260 (828.4 pb^-1), " \
  "703_4360 (543.9 pb^-1), 703_4420 (1043.9 pb^-1), 703_4600 (586.9 pb^-1).")

alg.apply(event_selection)

# --- Datasets: 5 main high-luminosity points ---
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")

data_points = [data_4230, data_4260, data_4360, data_4420, data_4600]

inc_mc_4230 = DatasetManager.inclusive_mc.find("703_4230")
inc_mc_4260 = DatasetManager.inclusive_mc.find("703_4260")
inc_mc_4360 = DatasetManager.inclusive_mc.find("703_4360")
inc_mc_4420 = DatasetManager.inclusive_mc.find("703_4420")
inc_mc_4600 = DatasetManager.inclusive_mc.find("703_4600")

inc_mc_points = [inc_mc_4230, inc_mc_4260, inc_mc_4360, inc_mc_4420, inc_mc_4600]

# Signal MC for scan
sig_mcs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_piD0Dstar"
  config.events        = 500_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg.execute_on(data_points + inc_mc_points + sig_mcs)