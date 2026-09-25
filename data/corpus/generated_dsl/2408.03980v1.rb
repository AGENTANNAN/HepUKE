# ============================================================================
# Datasets: D_s scan points (4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219,
# 4.226 GeV) plus the corresponding inclusive MC samples.
# Sample-name convention: [BOSS_version]_[CMS_energy_in_MeV]
# ============================================================================
data_4128 = DatasetManager.real_data.find("705_4130")   # 4.128 GeV
data_4157 = DatasetManager.real_data.find("705_4160")   # 4.157 GeV
data_4178 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4189 = DatasetManager.real_data.find("703_4190")   # 4.189 GeV
data_4199 = DatasetManager.real_data.find("703_4200")   # 4.199 GeV
data_4209 = DatasetManager.real_data.find("703_4210")   # 4.209 GeV
data_4219 = DatasetManager.real_data.find("703_4220")   # 4.219 GeV
data_4226 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV

incMC_4128 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4157 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4189 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4199 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4209 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4219 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

data_points = [data_4128, data_4157, data_4178, data_4189,
               data_4199, data_4209, data_4219, data_4226]
incMCs      = [incMC_4128, incMC_4157, incMC_4178, incMC_4189,
               incMC_4199, incMC_4209, incMC_4219, incMC_4226]

# ============================================================================
# Decay card: signal D_s+ -> gamma rho(770)+, rho+ -> pi+ pi0, pi0 -> gamma gamma,
# in e+e- -> D_s*+ D_s- events (D_s*+ -> gamma D_s+).
# ============================================================================
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s+
    1.000 gamma rho+ PHSP;
    Enddecay

    Decay rho+
    1.000 pi+ pi0 VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC, one sample per scan point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Ds_gamma_rho_exclusive_mc"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# ============================================================================
# Tag-based analysis (D_s tag): ST tag on D_s- + signal side D_s+ -> gamma rho+
# ============================================================================
alg_name = "DsGammaRho"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.178]})   # nominal scan energy
   .with_decay_card(decay_card_signal)

# --- Tag side: single tag on D_s- in five hadronic modes -------------------
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKKPiPi0, :DstoKsKPiPi, :DstoKPiPi
  t.charm(-1)                                   # tag the D_s- (ST)
end

# --- Signal side: D_s+ -> gamma rho+ -> gamma pi+ pi0 -> gamma pi+ gamma gamma
alg.signal_side do |s|
  s.photons 3..48              # 3-48 photons (event-level extra-shower cap)
  s.min_photon_energy 0.025    # >= 25 MeV
  s.min_photon_angle 10.0      # >= 10 degrees to charged tracks
  s.charged(pip: 1)            # exactly one pi+
  s.require_charge(1)          # net charge +1
end

# --- Kinematic fit: 4-momentum constraint + pi0 mass constraint ------------
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# --- BOSS-side selection with no dedicated DSL construct -------------------
alg.note(:recoil_mass_window,
  "energy-dependent window on the recoil mass against the tagged D_s- " \
  "(M_recoil ~ M(D_s*+)) used to select D_s*+ D_s- events; the window " \
  "boundaries depend on the CMS energy of each scan point")

alg.apply
alg.execute_on(data_points + incMCs + exMCs_signal)