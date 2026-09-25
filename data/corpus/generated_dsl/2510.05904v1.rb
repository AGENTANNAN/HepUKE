### Dataset description ###
data_4180 = DatasetManager.real_data.find("703_4180")       # primary real dataset at √s = 4178 MeV (4.128–4.226 GeV range)
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")   # corresponding inclusive MC sample

# Decay card for the signal process (Ds+ -> K0 mu+ nu_mu, K0 -> KS0 -> pi+ pi-)
decay_card_signal = <<~DECAYCARD
    Decay D_s+
    1.0000 K0 mu+ nu_mu PHSP;
    Enddecay

    Decay K0
    1.0000 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### Event selection (BOSS) — tag-based double-tag analysis ###
alg_name = "DsTagK0MuNu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.178]})   # √s = 4178 MeV
   # Inexpressible BOSS-side background suppressions (preserved for the systematic-uncertainties skill)
   .note(:background_veto_k0mu,
         "M(K0 mu+) < 1.70 GeV/c^2 required to suppress the Ds+ -> K0 pi+ background")
   .note(:background_veto_unused_photon,
         "maximum unused photon energy below 0.15 GeV required to suppress pi0 backgrounds")
   .note(:background_veto_dsstar,
         "Ds*± Ds∓ pair constraint (chi2 < 40) applied to suppress combinatorial background")
   .note(:pid_correction_method,
         "muon identification uses the likelihood method with L_mu > 0.001, L_mu > L_e, L_mu > L_K and an EMC energy deposit 0.1 < E < 0.3 GeV; this deviates from the fixed signal-side lepton PID recipe and must be applied in the generated code")
   .with_decay_card(decay_card_signal)

# Tag side: Ds- reconstructed through the 14 hadronic modes (one tag candidate per event = single tag)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,          # K- K+ pi-
          :DstoKKPiPi0,       # K- K+ pi- pi0
          :DstoKsK,           # KS0 K-
          :DstoKsKPi0,        # KS0 K- pi0
          :DstoKsPi,          # KS0 pi-
          :DstoKsPiPi0,       # KS0 pi- pi0
          :DstoKPiPi,         # K- pi+ pi-
          :DstoKPiPiPi0,      # K- pi+ pi- pi0
          :DstoPiPiPi,        # pi- pi+ pi-
          :DstoPiPiPiPi0,     # pi- pi+ pi- pi0
          :DstoPiPiPiPiPi,    # pi- pi+ pi- pi+ pi-
          :DstoPiEta,         # pi- eta
          :DstoPiPi0,         # pi- pi0
          :DstoPiEtaPiPiPi0   # pi- eta pi+ pi- pi0
  t.charm(-1)   # pin the tagged side to Ds-
end

# Signal side: one pi+, one pi- (from KS0 -> pi+ pi-) and one mu+; one missing nu_mu
alg.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)   # mu+ via the lepton PID recipe
  s.require_charge 1                  # total visible signal-side charge (Ds+)
  s.missing :nu_mu                    # massless missing neutrino (semileptonic tag), mass fixed to zero
end

# Kinematic fit: constrain the Ds four-momentum; loose chi2 cut at BOSS level
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200   # paper's chi2 < 40 is applied later at ROOT level
end

alg.apply
root_files = alg.execute_on([data_4180, incMC_4180])