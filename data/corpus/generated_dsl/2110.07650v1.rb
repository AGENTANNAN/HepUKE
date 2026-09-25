### Dataset description ###
# Real data at the six energy points (BOSS 7.0.3 samples, named <BOSS>_<Ecms(MeV)>)
data_4178 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4189 = DatasetManager.real_data.find("703_4190")   # 4.189 GeV
data_4199 = DatasetManager.real_data.find("703_4200")   # 4.199 GeV
data_4209 = DatasetManager.real_data.find("703_4210")   # 4.209 GeV
data_4219 = DatasetManager.real_data.find("703_4220")   # 4.219 GeV
data_4226 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV

# Inclusive MC at the same six energy points
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4189 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4199 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4209 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4219 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

datasets = [data_4178, data_4189, data_4199, data_4209, data_4219, data_4226,
            incMC_4178, incMC_4189, incMC_4199, incMC_4209, incMC_4219, incMC_4226]

### Decay card for the signal process (EvtGen format, EvtGen particle names) ###
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000  D_s+  D_s-        PHSP;
    Enddecay

    Decay D_s+
    1.000  K_S0  K_S0  pi+   PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-          PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-       PHSP;
    Enddecay

    End
DECAYCARD

### Tag-analysis specification ###
# Double-tag technique: tag-side Ds- reconstructed from pre-stored DTag candidates,
# signal-side Ds+ -> K_S0 K_S0 pi+ reconstructed from the tracks the tag did not use.
# A Selection object is NOT used here; the tag blocks are the selection.
alg_name = "DsToKsKsPiAmp"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })
   .with_decay_card(decay_card_signal)
   # 8C fit = 4C + 2xK_S0 mass + tag Ds mass + signal Ds mass + Ds* mass.
   # The generated code imposes the 4C and Ds-mass constraints; the two K_S0 and
   # the Ds* mass constraints are added by the BOSS job options.
   .note(:kinematic_fit_constraints,
     "8C kinematic fit: 4C (four-momentum conservation) + 2x K_S0 mass + tag Ds mass + " \
     "signal Ds mass + Ds* mass. Generated code applies the 4C and Ds-mass constraints; " \
     "the two K_S0 mass constraints and the Ds*-mass constraint are supplied in the BOSS job options.")
   .note(:signal_mass_window,
     "signal Ds+ invariant-mass window [1.950, 1.990] GeV/c^2 applied on the reconstructed " \
     "Ds+ -> K_S0 K_S0 pi+ candidate")
   .note(:background_veto,
     "the branching-fraction-measurement double-tag selection (no kinematic fit) requires the " \
     "soft pion momentum > 0.1 GeV/c to suppress D*+ decays")
   .note(:eta_prime_proxy,
     "the pi- eta' tag mode is handled through a pi- eta proxy")

# Tag side: Ds- reconstructed in the eight hadronic modes.
# The two charge-conjugate K_S0 K pi pi states (K_S0 K- pi- pi+ and K_S0 K+ pi- pi-) share
# one decay symbol.
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,        # K_S0 K-
          :DstoKKPi,       # K+ K- pi-
          :DstoKKPiPi0,    # K+ K- pi- pi0
          :DstoKsKPiPi,    # K_S0 K- pi- pi+ and its charge conjugate
          :DstoPiPiPi,     # pi- pi- pi+
          :DstoPiEtaPrime, # pi- eta'
          :DstoKPiPi       # K- pi+ pi-
  t.charm -1               # pin the tagged side to Ds-
end

# Signal side: Ds+ -> K_S0 K_S0 pi+ ; K_S0 -> pi+ pi-, so the final charged tracks are 3 pi+ + 2 pi-
alg.signal_side do |s|
  s.charged(pip: 3, pim: 2)
end

# Kinematic fit over the derived participants
alg.fit do |f|
  f.constrain_four_momentum                                      # 4C energy-momentum conservation
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D_s)  # tag Ds mass constraint
  f.chi2_cut 200
end

# Reconstruct the tag locally; beam energy is read from the conditions DB
alg.dtag_reconstruction do |d|
  d.local true
  d.beam_energy :db
end

alg.apply
alg.execute_on(datasets)