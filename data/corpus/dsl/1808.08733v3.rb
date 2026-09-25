# 1808.08733v3: Precision measurements of e+e- → Ks0 K± π∓ Born cross sections
# at center-of-mass energies between 3.8 and 4.6 GeV at BESIII
# Multi-energy ConExc scan with ISR correction (Born cross section measurement)

# === Datasets: BOSS 703 scan points matching the paper's 15 energy points ===
d_3810 = DatasetManager.real_data.find("703_3810")    # √s ≈ 3.808 GeV
d_3900 = DatasetManager.real_data.find("703_3900")    # √s ≈ 3.896 GeV
d_4009 = DatasetManager.real_data.find("703_4009")    # √s ≈ 4.009 GeV
d_4090 = DatasetManager.real_data.find("703_4090")    # √s ≈ 4.086 GeV
d_4190 = DatasetManager.real_data.find("703_4190")    # √s ≈ 4.189 GeV
d_4200 = DatasetManager.real_data.find("703_4200")    # √s ≈ 4.200 GeV
d_4210 = DatasetManager.real_data.find("703_4210")    # √s ≈ 4.208 GeV
d_4230 = DatasetManager.real_data.find("703_4230")    # √s ≈ 4.226 GeV
d_4245 = DatasetManager.real_data.find("703_4245")    # √s ≈ 4.242 GeV
d_4260 = DatasetManager.real_data.find("703_4260")    # √s ≈ 4.258 GeV
d_4310 = DatasetManager.real_data.find("703_4310")    # √s ≈ 4.308 GeV
d_4360 = DatasetManager.real_data.find("703_4360")    # √s ≈ 4.358 GeV
d_4390 = DatasetManager.real_data.find("703_4390")    # √s ≈ 4.387 GeV
d_4420 = DatasetManager.real_data.find("703_4420")    # √s ≈ 4.416 GeV
d_4600 = DatasetManager.real_data.find("703_4600")    # √s ≈ 4.600 GeV

scan_datasets = [d_3810, d_3900, d_4009, d_4090, d_4190, d_4200, d_4210, d_4230,
                 d_4245, d_4260, d_4310, d_4360, d_4390, d_4420, d_4600]

# === ConExc decay card for e+e- → K_S0 K± π∓ (ISR Born cross section) ===
# ConExc mode 9 = K_S K+ π-, mode 10 = K_S K- π+ (charge-conjugate)
# Equal BR for both modes; vhdr final state used for reconstruction particle list
decay_card = <<~DECAYCARD
  Decay vpho
  0.5 ConExc 9;
  0.5 ConExc 10;
  Enddecay
  Decay vhdr
  1 K_S0 K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# === Exclusive signal MC for the full energy scan ===
sig_mc = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name   = "sig_Ks0Kpi_conexc"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# === Event selection (BOSS) ===
alg = Algorithm.new("Ks0KpiCrossSection")
alg.set_header(["Ks0KpiCrossSectionAlg/Ks0KpiCrossSection.h"])

sel = Selection.new

# Charged track selection: 2π+2π- (Ks0→π+π- + K±π∓), net charge zero
# |cosθ| < 0.93, |Vz| < 10 cm, |Vr| < 1 cm
sel.select_track do
  cos_theta   0.93
  Vz          100.0
  Vr          10.0
  nChrp       "==2"
  nChrn       "==2"
  nNet        "==0"
end

# Photon selection (no photons required; select for EMC-aware tracking)
sel.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
end

# Particle identification using dE/dx + TOF probability method
# Kaon and pion hypotheses with 0.001 probability cut
sel.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
end

# Suppress photon conversion background: Ks0 pions must satisfy E/p < 0.8
sel.for_each(:pip) { where { eraw / p < 0.8 }; remove }
sel.for_each(:pim) { where { eraw / p < 0.8 }; remove }

# Reconstruct Ks0 → π+π- with secondary vertex fit
# Decay length > 2σ enforced by vertex fit; best candidate by mass difference
sel.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# 4C kinematic fit under e+e- → K_S0 K± π∓ hypothesis
# Loose chi2 cut (200); optimal tight cut determined in ROOT
sel.kinematic_fit([:K_S0, :kp, :km, :pip, :pim]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg.with_decay_card(decay_card)
   .note(:ks0_mass_window,
     "|M(π+π-) - M(K_S0)| < 0.020 GeV/c²; signal yield via event counting " \
     "in signal region with sideband subtraction in ROOT")
   .note(:ks0_sideband,
     "Sideband regions: m(π+π-) ∈ (0.435, 0.455) ∪ (0.545, 0.565) GeV/c²")
   .note(:ks0_vertex_cuts,
     "Ks0 pions: distance of closest approach < 25 cm (z) and 20 cm (rφ); " \
     "decay length > 2σ of vertex fit; multiple Ks0 candidates resolved by smallest χ²")
   .note(:conexc_isr,
     "ISR correction factor (1+δ_ISR) by iterative procedure converging to 1.0%; " \
     "BABAR Born cross sections used as initial input")
   .note(:charge_conjugate,
     "Both K_S0 K+ π- and charge-conjugate K_S0 K- π+ included; " \
     "ConExc modes 9+10 with equal branching fraction")
   .note(:pwa_efficiency,
     "Detection efficiency from PHSP MC reweighted by partial-wave analysis; " \
     "intermediate states: K*(892), K2*(1430), K3*(1780), a2(1320), ρ(1700), ρ(2150)")
   .note(:systematic_uncertainties,
     "Tracking 2.0%, PID 2.0%, Ks0 reconstruction 1.2%, kinematic fit 0.5%, " \
     "signal model 2.0%, signal yield 1.8%, ISR factor 1.0%, luminosity 1.0%, BF 0.1%; " \
     "total 4.4%")
   .apply(sel)

alg.execute_on(scan_datasets + sig_mc)