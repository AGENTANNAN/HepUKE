### Dataset description ###
# Six c.m. energy points: 4178, 4189, 4199, 4209, 4219, 4226 MeV
data_4178 = DatasetManager.real_data.find("703_4180")     # e+e- data at 4178 MeV
data_4189 = DatasetManager.real_data.find("703_4190")     # 4189 MeV
data_4199 = DatasetManager.real_data.find("703_4200")     # 4199 MeV
data_4209 = DatasetManager.real_data.find("703_4210")     # 4209 MeV
data_4219 = DatasetManager.real_data.find("703_4220")     # 4219 MeV
data_4226 = DatasetManager.real_data.find("703_4230")     # 4226 MeV

incMC_4178 = DatasetManager.inclusive_mc.find("703_4180") # inclusive MC at the same energies
incMC_4189 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4199 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4209 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4219 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

### Decay card for the signal process (EvtGen format) ###
# e+e- -> D_s+ D_s- followed by D_s+ -> tau+ nu_tau, tau+ -> e+ nu_e nu_tau.
# The D_s- (tag side) is left to the generator's default decay table so the
# 11 hadronic tag modes are populated inclusively.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s- PHSP;
    Enddecay

    Decay D_s+
    1.0000 tau+ nu_tau PHSP;
    Enddecay

    Decay tau+
    1.0000 e+ nu_e nu_tau PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC (same card, one sample per energy point) ###
data_points = [data_4178, data_4189, data_4199, data_4209, data_4219, data_4226]
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ds_taunu_signal_mc"   # auto-suffixed per energy point
  config.events        = 500000                 # 500k events per energy point
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based D_s analysis ###
alg = TagAnalysis.new("DsTauNu")
alg.set_header(["DsTauNuAlg/DsTauNu.h"])
   .set_constant({"ECMS" => [:double, 4.178]})   # nominal; per-run energy read from DB

# Tag side: reconstruct the D_s- in 11 hadronic modes (single tag)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,                          # K_S K-
          :DstoKKPi,                          # K+ K- pi-
          :DstoKKPiPi0,                       # K+ K- pi- pi0
          :DstoKsKmPiPi,                      # K_S K- pi+ pi-
          :DstoKsKpPiPi,                      # K_S K+ pi- pi-
          :DstoPiPiPi,                        # pi+ pi- pi-
          :DstoPiEta,                         # pi- eta
          :DstoPiPi0Eta,                      # pi- pi0 eta
          :DstoPiEtaPrimeToPiPiEta,           # pi- eta'(-> pi+ pi- eta)
          :DstoPiEtaPrimeToGammaRho,          # pi- eta'(-> gamma rho0)
          :DstoKPiPi                          # K- pi+ pi-
  t.charm -1                                  # pin the tagged side to D_s-
  t.window :mBC, min: 1.89, max: 2.04         # single-tag M_ST signal region
end

# Signal side: what the tag did not use
alg.signal_side do |s|
  s.charged(ep: 1)        # exactly one extra positron
  s.require_charge(1)     # charge opposite to the D_s- tag
  s.missing :nu_tau       # massless missing nu_tau
  s.min_photon_angle 5.0  # accept nearby showers to recover FSR photons
end

# 4C kinematic fit: tag + e+ + nu_tau constrained to the lab four-momentum
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures that cannot be expressed in formal DSL syntax
alg.note(:recoil_mass_window, "energy-point-dependent recoil-mass windows applied to the tag for D_s*D_s / D_sD_s topologies; e.g. 2.050-2.195 GeV at 4178 MeV and 2.040-2.220 GeV at 4226 MeV; stored per energy point and windowed in ROOT")
   .note(:pid_correction_method, "signal positron selected with CL_e > 0.1%, CL_e/(CL_e+CL_pi+CL_K) > 0.8, p > 0.2 GeV/c, and E/p > 0.8")
   .note(:fsr_recovery, "photons within 5 degrees of the signal positron are added to recover FSR")
   .note(:background_veto, "Bhabha events suppressed by a relative-probability sum cut (<2.0 for n>1, <0.9 for n=1)")
   .note(:background_sources, "backgrounds from non-D_s-, D_s+ -> K_L0 e+ nu_e, and D_s+ -> X e+ nu_e accounted for")
   .note(:extra_energy_cut, "double-tag yield extracted from E_extra^tot < 0.4 GeV")

# Render the tag specification (apply takes no Selection argument)
alg.with_decay_card(decay_card_signal)
alg.apply

# Execute on real data, inclusive MC, and exclusive signal MC
alg.execute_on(data_points +
               [incMC_4178, incMC_4189, incMC_4199, incMC_4209, incMC_4219, incMC_4226] +
               exMCs_signal)