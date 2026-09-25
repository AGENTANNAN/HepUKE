# =============================================================================
# Ds+ -> mu+ nu_mu  --  tag-based (single-tag) measurement
# 8 continuum energy points, 4.128 - 4.226 GeV  (7.33 fb^-1)
# Same selection chain is applied to data, inclusive MC and exclusive signal MC
# =============================================================================

### Dataset preparation ###

# Real data: eight continuum points
data_points = [
  DatasetManager.real_data.find("705_4130"),  # 4.128 GeV
  DatasetManager.real_data.find("705_4160"),  # 4.157 GeV
  DatasetManager.real_data.find("703_4180"),  # 4.178 GeV
  DatasetManager.real_data.find("703_4190"),  # 4.189 GeV
  DatasetManager.real_data.find("703_4200"),  # 4.199 GeV
  DatasetManager.real_data.find("703_4210"),  # 4.209 GeV
  DatasetManager.real_data.find("703_4220"),  # 4.219 GeV
  DatasetManager.real_data.find("703_4230")   # 4.226 GeV
]

# Corresponding inclusive MC at each energy point
inc_mc_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Decay card for the exclusive signal MC:
#   e+e- -> Ds*+ Ds-  (continuum top-mother convention: psi(4260) + KKMC)
#   Ds+ -> mu+ nu_mu        (signal side)
#   Ds- -> K+ K- pi-        (tag side)
#   Ds* -> Ds gamma         (transition photon)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000  D_s*+  D_s-   PHSP;
    Enddecay

    Decay D_s*+
    1.000  D_s+  gamma   PHSP;
    Enddecay

    Decay D_s+
    1.000  mu+  nu_mu    PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-   PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal sample, one per energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ds_munu_tag_exmc"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###

alg_name = "DsTagMuNu"
alg = TagAnalysis.new(alg_name)
# ECMS is only a reference value; the tag fit reads the per-run beam energy from the DB.
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.180] })
   .with_decay_card(decay_card_signal)

# --- Tag side: Ds- reconstructed in the hadronic modes (single tag) ---
#   K+K-pi- , K+K-pi-pi0 , pi+pi-pi- , KS K- , KS K-pi0 , KS KS pi- , KS K-pi+pi- ,
#   eta(gg)pi- , eta(3pi)pi- , eta'(gg)pi+pi-pi- , eta'(g pi+pi-)pi+pi-pi- ,
#   eta(gg)rho , eta(3pi)rho
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi,
          :DstoKsK, :DstoKsKPi0, :DstoKsKsPi, :DstoKsKPiPi,
          :DstoEtaPiEtaToGG,
          :DstoEtaPiEtaToPiPiPi0,
          :DstoEtaPrimePiPiPiEtaPrimeToGG,
          :DstoEtaPrimePiPiPiEtaPrimeToRhoGam,
          :DstoEtaRhoEtaToGG,
          :DstoEtaRhoEtaToPiPiPi0
end

# --- Signal side: mu+ + transition photon + undetected nu_mu ---
alg.signal_side do |s|
  s.charged(mup: 1)          # exactly one signal-side mu+ candidate
  s.photons 1                # at least one photon from the Ds* -> Ds gamma transition
  s.min_photon_energy 0.025  # photon energy above 25 MeV
  s.min_photon_angle 10.0    # photon isolated from charged tracks
  s.missing :nu_mu           # the nu_mu is undetected (massless missing particle)
end

# --- 4C kinematic fit: four-momentum conservation + Ds and Ds* mass constraints ---
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)                       # Ds mass (tag side)
  f.invariant_mass_of(:mup, :nu_mu, :gamma).constrain_to_nominal_mass_of(:Ds_star)   # Ds* mass (signal side)
  f.chi2_cut 200
end

# --- BOSS-side procedures that cannot be expressed in formal DSL syntax ---
alg.note(:tag_modes_unavailable,
         "two additional Ds- tag modes of the analysis are absent from the DTagAlg
          tag-mode vocabulary and are therefore not declared in tag_side")
alg.note(:pid_correction_method,
         "signal mu+ candidate required to have EMC deposited energy between 0 and 0.3 GeV
          together with muon-chamber hit-depth conditions; not expressible by the
          signal-side charged(mup: 1) declaration alone")
alg.note(:photon_selection_method,
         "the Ds* -> Ds gamma transition photon is chosen by the minimum |deltaE| method")
alg.note(:background_veto,
         "the maximum energy of extra (non-signal) photons is required to be below 0.3 GeV")

alg.apply
root_files = alg.execute_on(data_points + inc_mc_points + exMCs_signal)