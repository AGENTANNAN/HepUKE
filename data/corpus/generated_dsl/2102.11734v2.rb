# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Six centre-of-mass energies: 4.178, 4.190, 4.200, 4.210, 4.220 and 4.226 GeV (6.32 fb^-1 total)
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4.190 GeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4.200 GeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4.210 GeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4.220 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

data_points  = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]
incMC_points = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

# Decay card for the muonic signal: e+e- -> Ds+ Ds- at psi(4260); signal Ds+ -> mu+ nu_mu.
# The tag-side Ds- is left undecayed in the card (EvtGen default table -> inclusive decay).
decay_card_munu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.0000 mu+ nu_mu PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for the tauonic signal: Ds+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau;
# tag-side Ds- again left inclusive.
decay_card_taunu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.0000 tau+ nu_tau PHSP;
  Enddecay

  Decay tau+
  1.0000 pi+ anti-nu_tau PHSP;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive signal MC per leptonic mode, one sample per energy point
exMC_munu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_munu"
  config.events        = 500_000
  config.decay_card    = decay_card_munu
  config.cross_section = :default
end

exMC_taunu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_taunu"
  config.events        = 500_000
  config.decay_card    = decay_card_taunu
  config.cross_section = :default
end

### Muonic channel: Ds+ -> mu+ nu_mu (tag Ds- in 13 hadronic modes) ###
alg_munu = TagAnalysis.new("DsToMuNu")
alg_munu.set_header(["DsToMuNuAlg/DsToMuNu.h"])
        .with_decay_card(decay_card_munu)
        .note(:tag_track_selection,
              "DTagAlg tag-side charged tracks: |cos(theta)|<0.93, |Vz|<10 cm, |Vxy|<1 cm; " \
              "K_S0 and Lambda daughter tracks are exempt from these cuts")
        .note(:tag_photon_selection,
              "DTagAlg tag-side photons: EMC clusters with E>25 MeV in the barrel (|cos(theta)|<0.80) " \
              "or E>50 MeV in the endcap (0.86<|cos(theta)|<0.92); EMC time in (0,700) ns; " \
              ">10 degrees from any charged track")
        .note(:tag_pi0_reconstruction,
              "pi0 -> gamma gamma pairs with 0.115<M<0.150 GeV passing a 1C mass-constrained fit chi2<20")
        .note(:tag_ks0_reconstruction,
              "K_S0 -> pi+ pi- pairs with a secondary-vertex fit, positive decay length, " \
              "0.485<M<0.510 GeV")
        .note(:pid_correction_method,
              "tag-side K/pi separation from combined dE/dx and TOF confidence levels: " \
              "kaon if CL_K>CL_pi, pion if CL_pi>CL_K")
        .note(:tag_deltaE_window,
              "tag DeltaE cut is mode-dependent (about +/-3 sigma around zero) and cannot be " \
              "expressed as a single DSL window value; applied as a stored-variable window in ROOT")
        .note(:signal_muon_pid,
              "signal muon PID: L'_mu>0.001, L'_mu>L'_e and L'_mu>L'_K")
        .note(:signal_extra_energy,
              "no extra photon beyond the Ds* radiative photon: E_extra<0.3 GeV")
        .note(:dsstar_photon,
              "Ds*+ -> gamma Ds+ radiative photon selected in the Ds* rest frame with energy " \
              "in 119-149 MeV")

# Tag side: Ds- reconstructed in 13 hadronic tag modes (single tag; the signal Ds+ is
# inferred from the recoil / missing mass). Explicit tag m_BC > 2.05 GeV window requested.
alg_munu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,        # Ds- -> K+ K- pi-
          :DstoKPiPi,       # Ds- -> K+ pi- pi-
          :DstoPiPiPi,      # Ds- -> pi+ pi- pi-
          :DstoKKPiPi0,     # Ds- -> K+ K- pi- pi0
          :DstoKPiPiPi0,    # Ds- -> K+ pi- pi- pi0
          :DstoPiPiPiPi0,   # Ds- -> pi+ pi- pi- pi0
          :DstoKKPi0,       # Ds- -> K+ K- pi0
          :DstoKPiPi0,      # Ds- -> K+ pi- pi0
          :DstoPiPiPi0,     # Ds- -> pi+ pi- pi0
          :DstoKKPiPiPi,    # Ds- -> K+ K- pi- pi+ pi-
          :DstoKKPiPi0Pi0,  # Ds- -> K+ K- pi- pi0 pi0
          :DstoKPiPiPiPi,   # Ds- -> K+ pi- pi- pi+ pi-
          :DstoPiPiPiPiPi   # Ds- -> pi+ pi- pi- pi+ pi-
  t.charm -1                # pin the reconstructed (tagged) side to Ds-
  t.window :mBC, min: 2.05  # M_BC > 2.05 GeV/c^2
end

# Signal side: exactly one mu+ opposite the tag, the Ds* radiative photon, and the missing nu_mu.
alg_munu.signal_side do |s|
  s.charged(mup: 1)         # exactly one signal muon
  s.photons 1               # the Ds*+ -> gamma Ds+ radiative photon
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_mu          # missing (massless) neutrino
end

# 4-momentum conservation including the missing neutrino.
alg_munu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_munu.apply   # no Selection argument for a TagAnalysis
root_files_munu = alg_munu.execute_on(data_points + incMC_points + exMC_munu)

### Tauonic channel: Ds+ -> tau+ nu_tau with tau+ -> pi+ anti-nu_tau ###
alg_taunu = TagAnalysis.new("DsToTauNu")
alg_taunu.set_header(["DsToTauNuAlg/DsToTauNu.h"])
         .with_decay_card(decay_card_taunu)
         .note(:tag_track_selection,
               "DTagAlg tag-side charged tracks: |cos(theta)|<0.93, |Vz|<10 cm, |Vxy|<1 cm; " \
               "K_S0 and Lambda daughter tracks are exempt from these cuts")
         .note(:tag_photon_selection,
               "DTagAlg tag-side photons: EMC clusters with E>25 MeV in the barrel (|cos(theta)|<0.80) " \
               "or E>50 MeV in the endcap (0.86<|cos(theta)|<0.92); EMC time in (0,700) ns; " \
               ">10 degrees from any charged track")
         .note(:tag_pi0_reconstruction,
               "pi0 -> gamma gamma pairs with 0.115<M<0.150 GeV passing a 1C mass-constrained fit chi2<20")
         .note(:tag_ks0_reconstruction,
               "K_S0 -> pi+ pi- pairs with a secondary-vertex fit, positive decay length, " \
               "0.485<M<0.510 GeV")
         .note(:pid_correction_method,
               "tag-side K/pi separation from combined dE/dx and TOF confidence levels: " \
               "kaon if CL_K>CL_pi, pion if CL_pi>CL_K")
         .note(:tag_deltaE_window,
               "tag DeltaE cut is mode-dependent (about +/-3 sigma around zero) and cannot be " \
               "expressed as a single DSL window value; applied as a stored-variable window in ROOT")
         .note(:signal_pion_pid,
               "signal pion: CL_pi>CL_K and CL_pi>CL_p, with muon veto L'_mu<L'_pi")
         .note(:signal_extra_energy,
               "no extra photon beyond the Ds* radiative photon: E_extra<0.3 GeV")
         .note(:dsstar_photon,
               "Ds*+ -> gamma Ds+ radiative photon selected in the Ds* rest frame with energy " \
               "in 119-149 MeV")

alg_taunu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKPiPi,
          :DstoPiPiPi,
          :DstoKKPiPi0,
          :DstoKPiPiPi0,
          :DstoPiPiPiPi0,
          :DstoKKPi0,
          :DstoKPiPi0,
          :DstoPiPiPi0,
          :DstoKKPiPiPi,
          :DstoKKPiPi0Pi0,
          :DstoKPiPiPiPi,
          :DstoPiPiPiPiPi
  t.charm -1
  t.window :mBC, min: 2.05
end

# Signal side: exactly one pi+ opposite the tag, the Ds* radiative photon, and the missing neutrino(s).
alg_taunu.signal_side do |s|
  s.charged(pip: 1)         # exactly one signal pion
  s.photons 1               # the Ds*+ -> gamma Ds+ radiative photon
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu             # combined missing (massless) neutrino system
end

# 4-momentum conservation including the missing neutrino(s).
alg_taunu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_taunu.apply  # no Selection argument for a TagAnalysis
root_files_taunu = alg_taunu.execute_on(data_points + incMC_points + exMC_taunu)