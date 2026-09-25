# ============================================================
# Dataset preparation
# ============================================================
# Six BESIII c.m. energy points from 4.178 to 4.226 GeV (4180-4230 MeV), BOSS 703
data_points = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230")
]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Decay card for the signal process: psi(4260) -> gamma D_s+ D_s-
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 K_S0 K- pi+ pi+ PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# 1M-event exclusive signal MC, generated for each energy point of the scan
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_psip4260_gDsDs"
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

# ============================================================
# Tag analysis: single D_s- tag, signal side D_s*+ -> gamma D_s+
# ============================================================
alg = TagAnalysis.new("DsTagDsStarGamma")
alg.set_header(["DsTagDsStarGammaAlg/DsTagDsStarGamma.h"])
   .set_constant({ "ECMS" => [:double, 4.226] })
   .set_alias({ "std::vector<double>" => "Vdouble" })
   # Helix-parameter correction applied to all charged tracks before the 4C fit
   .note(:helix_correction, "helix-parameter track correction applied to all charged
     tracks before the 4C kinematic fit; residual efficiency difference estimated by
     re-running the BOSS selection with and without the correction on signal MC")
   # D0 / anti-D0 mass veto
   .note(:background_veto, "events whose invariant mass is consistent with a
     D0 / anti-D0 candidate are vetoed to suppress D0(Dbar0) backgrounds; veto window
     taken from the D0 mass peak in inclusive MC")
   # Only eight of the nine D_s- hadronic tag modes are used
   .note(:tag_mode_availability, "eight of the nine D_s- hadronic tag modes are used;
     D_s- -> pi- eta' is unavailable and is therefore omitted")
   .with_decay_card(decay_card_signal)

# --- Tag side: single D_s- tag in eight hadronic modes ---
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,        # K_S0 K-
          :DstoKKPi,       # K+ K- pi-
          :DstoKsKPiPi,    # K_S0 K- pi+ pi-
          :DstoKPiPi,      # K+ pi- pi-
          :DstoPiEta,      # pi- eta
          :DstoKsKPi0,     # K_S0 K- pi0
          :DstoKKPiPiPi,   # K+ K- pi+ pi- pi-
          :DstoPiPiPi      # pi+ pi- pi-
  t.charm -1               # tag the D_s- side
end

# --- Signal side: D_s*+ -> gamma D_s+ -> gamma K_S0 K- pi+ pi+ (K_S0 -> pi+ pi-) ---
alg.signal_side do |s|
  s.photons 1                      # at least one photon (from the D_s*+ -> gamma D_s+)
  s.min_photon_energy 0.025        # ... with energy above 25 MeV
  s.charged(km: 1, pip: 3, pim: 1) # K- and pi+ pi+ from the D_s+; pi+ pi- from K_S0
  s.require_charge(1)              # net charge +1
end

# --- 4C kinematic fit, chi2 < 200 ---
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Render and run on real data, inclusive MC and the exclusive signal MC
alg.apply
alg.execute_on(data_points + incMC_points + exMC_signal)