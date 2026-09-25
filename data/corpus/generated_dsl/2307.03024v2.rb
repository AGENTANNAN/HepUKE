# =============================================================================
# e+e- -> Ds*+ Ds-  with  Ds*+ -> gamma/pi0 Ds+  and the semileptonic signal
# Ds+ -> K+ K- mu+ nu_mu ; the opposite-side Ds- is reconstructed as a tag
# from the pre-stored DTag candidates (14 hadronic modes).  Tag-based analysis.
# =============================================================================

### Dataset preparation ###
# --- eight real-data points of the 4.128 - 4.226 GeV scan ---
data_4128 = DatasetManager.real_data.find("705_4130")   # sqrt(s) = 4128.5 MeV
data_4157 = DatasetManager.real_data.find("705_4160")   # sqrt(s) = 4157.4 MeV
data_4178 = DatasetManager.real_data.find("703_4180")   # sqrt(s) = 4178.0 MeV
data_4188 = DatasetManager.real_data.find("703_4190")   # sqrt(s) = 4188.8 MeV
data_4198 = DatasetManager.real_data.find("703_4200")   # sqrt(s) = 4198.9 MeV
data_4209 = DatasetManager.real_data.find("703_4210")   # sqrt(s) = 4209.2 MeV
data_4218 = DatasetManager.real_data.find("703_4220")   # sqrt(s) = 4218.7 MeV
data_4226 = DatasetManager.real_data.find("703_4230")   # sqrt(s) = 4226.3 MeV
data_points = [data_4128, data_4157, data_4178, data_4188,
               data_4198, data_4209, data_4218, data_4226]

# --- inclusive MC at 4.178 GeV (about 40x the data luminosity) ---
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")

# --- decay card for the semileptonic signal mode (phase space) ---
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Ds*+ Ds-  PHSP;
    Enddecay

    Decay Ds*+
    1.0000 gamma Ds+ PHSP;
    Enddecay

    Decay Ds+
    1.0000 K+ K- mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

# --- 200k phase-space signal MC events generated at each energy point ---
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ds_kkmunu_signal_mc"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "DsTagKKMuNu"
ds_tag = TagAnalysis.new(alg_name)
ds_tag.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({"ECMS" => [:double, 4.178]})            # reference (inclusive-MC) energy point
      .set_alias({"std::vector<double>" => "Vdouble"})
      .with_decay_card(decay_card_signal)

# ---- tag side: opposite-side Ds- in the 14 hadronic modes ----
ds_tag.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,        # K+ K- pi-
          :DstoKPiPi,       # K+ pi- pi-
          :DstoPiPiPi,      # pi+ pi- pi-
          :DstoKKPiPi0,     # K+ K- pi- pi0
          :DstoPi0Rho,      # pi0 rho
          :DstoPiPi0Eta,    # pi pi0 eta
          :DstoKsKmPiPi,    # K_S K- pi+ pi-
          :DstoKsKpPiPi,    # K_S K+ pi- pi-
          :DstoEtaPi,       # eta pi
          :DstoKsKsPi,      # K_S K_S pi-
          :DstoEtaPiPiPi,   # eta pi pi pi
          :DstoPi0PiPiEta,  # pi0 pi pi eta
          :DstoKsKPi0,      # K_S K- pi0
          :DstoKsK          # K_S K-
  t.charm -1                # tag the Ds- side (opposite to the Ds*+ signal side)
end

# ---- signal side: K+ K- mu+ (net charge +1) plus one missing nu_mu ----
ds_tag.signal_side do |s|
  s.charged(kp: 1, km: 1, mup: 1)  # exact charged multiset -> no unused charged track left;
                                   # mu+ is selected with the combined dE/dx-TOF-EMC likelihood
                                   # criteria L_mu > L_K, L_mu > L_e, L_mu > 0.001 (v1 thresholds)
  s.require_charge 1               # net charge of the signal side is +1
  s.missing :nu_mu                 # semileptonic neutrino (massless missing particle)
  s.min_photon_energy 0.2          # veto extra photons with E_gamma > 0.2 GeV
end

# ---- 4C kinematic fit: energy-momentum conservation, chi2 < 200 ----
ds_tag.fit do |f|
  f.constrain_four_momentum                               # tag + signal + missing = measured CMS 4-momentum
  f.invariant_mass_of(:kp, :km, :mup).between(0.0, 1.75)  # reject Ds+ -> K+ K- pi+ (pi mis-identified as mu)
  f.chi2_cut 200
end

# ---- BOSS-side procedures with no dedicated DSL construct ----
ds_tag.note(:background_veto,
        "Ds+ -> K+ K- pi+ background further rejected by M(K+ K- nu_mu) > 1.30 GeV/c^2 " \
        "(computed from the fitted missing four-momentum and windowed on the stored " \
        "missing-mass variable Umiss at ROOT level); events with any unused pi0 " \
        "candidate (not consumed by the tag or by the signal selection) are rejected")
      .note(:dsstar_transition,
        "the Ds*+ transition (gamma or pi0, both soft: E_gamma ~ 0.14 GeV) is assigned " \
        "among the unused photons / photon pairs by minimising |deltaE| when several " \
        "candidates survive; the Ds* mass constraint is imposed in the ROOT fit, not in " \
        "this BOSS selection (BOSS applies the 4C fit only)")

ds_tag.apply
root_files = ds_tag.execute_on(data_points + [incMC_4178] + exMCs_signal)